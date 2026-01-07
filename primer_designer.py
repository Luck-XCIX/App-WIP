# ============================================================================== 
# --- IMPORTS --- 
# ==============================================================================
# --- Standard Library Imports ---
import re
from itertools import product
from functools import lru_cache
# --- Third-party Imports ---
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# ============================================================================== 
# --- CORE BIOINFORMATICS UTILITY FUNCTIONS --- 
# ==============================================================================
def calculate_gc_content(sequence):
    """
    Calculates the GC content of a DNA sequence.
    """
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return round(gc_count / len(sequence) * 100, 2)

def extract_kmers(sequence, k_range):
    """
    Yields all k-mers of specified lengths from a sequence.
    """
    if not sequence:
        return
    for k in k_range:
        if 0 < k <= len(sequence):
            for i in range(len(sequence) - k + 1):
                yield sequence[i:i+k]

COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def are_bases_complementary(base1, base2):
    """Checks if two DNA bases are complementary using a fast lookup."""
    return COMPLEMENT_MAP.get(base1) == base2

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence using Biopython."""
    return str(Seq(sequence).reverse_complement())

def find_longest_complementary_run(seq1, seq2):
    """
    Finds the length of the longest perfectly complementary substring between two sequences.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return 0
    dp_row = [0] * (n + 1)
    max_len = 0
    for i in range(1, m + 1):
        prev_diagonal_val = 0
        for j in range(1, n + 1):
            current_dp_j = dp_row[j]
            if are_bases_complementary(seq1[i-1], seq2[j-1]):
                dp_row[j] = prev_diagonal_val + 1
                if dp_row[j] > max_len:
                    max_len = dp_row[j]
            else:
                dp_row[j] = 0
            prev_diagonal_val = current_dp_j
    return max_len

# ==============================================================================
# PERFORMANCE OPTIMIZATION: Global Tm cache
# ==============================================================================
_TM_CACHE = {}
_CACHE_HITS = 0
_CACHE_MISSES = 0

def calculate_melting_temp(sequence, na_conc=50e-3):
    """
    Calculates the melting temperature (Tm) using Biopython's nearest-neighbor model.
    NOW WITH CACHING: Reuses previously calculated Tm values for the same sequence.
    """
    global _CACHE_HITS, _CACHE_MISSES
    
    # Use uppercase for cache key consistency
    cache_key = (sequence.upper(), na_conc)
    
    if cache_key in _TM_CACHE:
        _CACHE_HITS += 1
        return _TM_CACHE[cache_key]
    
    _CACHE_MISSES += 1
    try:
        tm_value = round(mt.Tm_NN(str(sequence), Na=na_conc), 1)
        _TM_CACHE[cache_key] = tm_value
        return tm_value
    except Exception:
        return 0.0


def get_cache_stats():
    """Returns Tm cache statistics for debugging."""
    total = _CACHE_HITS + _CACHE_MISSES
    hit_rate = (_CACHE_HITS / total * 100) if total > 0 else 0
    return {
        "hits": _CACHE_HITS,
        "misses": _CACHE_MISSES,
        "total": total,
        "hit_rate": hit_rate,
        "cache_size": len(_TM_CACHE)
    }


def clear_tm_cache():
    """Clear the Tm cache (useful between different runs)."""
    global _TM_CACHE, _CACHE_HITS, _CACHE_MISSES
    _TM_CACHE.clear()
    _CACHE_HITS = 0
    _CACHE_MISSES = 0

def extract_isoform_id(name):
    """
    Extracts the isoform ID from a sequence name like 'CASC15-205 cDNA (4297 bp)'.
    Returns '205' from 'CASC15-205'.
    """
    match = re.search(r"CASC15-(\d+)", name)
    if not match:
        print(f"⚠️ Could not extract isoform ID from: {name}")
    return match.group(1) if match else ""

def check_primer_hits(primer, sequences):
    """
    Finds all sequence records that contain a given primer sequence
    or its reverse complement.
    """
    primer_rc = reverse_complement(primer)
    hits = [
        rec.id
        for rec in sequences
        if primer in str(rec.seq) or primer_rc in str(rec.seq)
    ]
    return {"total": len(hits), "isoforms": ",".join(hits)}

def find_candidate_primers(sequence, k_range, gc_min, gc_max, tm_min, tm_max):
    """
    Find candidate primers that meet GC and Tm criteria.
    OPTIMIZED: Early filtering - checks GC before calculating expensive Tm.
    """
    candidates = []
    gc_rejected = 0
    tm_rejected = 0
    
    for k in k_range:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            
            # OPTIMIZATION: Check GC first (cheap operation)
            gc = calculate_gc_content(kmer)
            if not (gc_min <= gc <= gc_max):
                gc_rejected += 1
                continue  # Skip expensive Tm calculation
            
            # Only calculate Tm if GC passed (expensive operation)
            tm = calculate_melting_temp(kmer)
            if tm_min <= tm <= tm_max:
                candidates.append(kmer)
            else:
                tm_rejected += 1
    
    # Debug info (can be removed in production)
    # print(f"  Filtered: {gc_rejected} by GC, {tm_rejected} by Tm | Accepted: {len(candidates)}")
    
    return candidates

# ============================================================================== 
# --- PRIMER DESIGN CORE LOGIC --- 
# ==============================================================================
@lru_cache(maxsize=10000)
def check_cached_primer_hits(primer, sequences_tuple):
    """
    A cached wrapper for `check_primer_hits` to avoid re-computing hits
    for the same primer. The `lru_cache` decorator memoizes the results.
    """
    sequences = [SeqIO.SeqRecord(Seq(seq), id=name) for name, seq in sequences_tuple]
    return check_primer_hits(primer, sequences)

def prepare_sequences_for_caching(sequences):
    """Converts a list of SeqRecord objects to a hashable tuple for caching."""
    return tuple((rec.id, str(rec.seq)) for rec in sequences)

def design_qpcr_primers(sequences, params, mode):
    """
    qPCR primer design scanning the full sequence.
    Specificity mode: Uses progressive relaxation (1 → MAX_RELAX hits) to find best primers.
    """
    MAX_RELAX = 4  # Maximum number of isoforms a primer can hit
    
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    
    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()
        
        if len(seq) < params["prod_min"]:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Find all candidate primers from full sequence
        fw_candidates = find_candidate_primers(
            seq, params["k_range"], params["gc_min"], params["gc_max"], 
            params["tm_min"], params["tm_max"]
        )
        
        rev_template_candidates = find_candidate_primers(
            seq, params["k_range"], params["gc_min"], params["gc_max"], 
            params["tm_min"], params["tm_max"]
        )
        
        # Map: actual reverse primer → template sequence
        rev_candidates_map = {
            reverse_complement(p).upper(): p.upper() 
            for p in rev_template_candidates
        }
        
        if not fw_candidates or not rev_candidates_map:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Check how many isoforms each primer hits
        fw_checks = {
            fw.upper(): check_cached_primer_hits(fw.upper(), sequences_tuple) 
            for fw in fw_candidates
        }
        rev_checks = {
            rv.upper(): check_cached_primer_hits(rv.upper(), sequences_tuple) 
            for rv in rev_candidates_map.keys()
        }
        
        # Convert to uppercase for consistency
        all_fw = [fw.upper() for fw in fw_candidates]
        all_rev = list(rev_candidates_map.keys())
        
        found_pair = None
        found_other_pairs = []
        
        if mode == "specificity":
            # PROGRESSIVE RELAXATION: Try 1, 2, 3, ... MAX_RELAX isoform hits
            for allowed_hits in range(1, MAX_RELAX + 1):
                print(f"  🔍 [{name}] Trying relaxation level: {allowed_hits} hit(s)...")
                
                # Filter primers that hit <= allowed_hits AND include target isoform
                isoform_id = extract_isoform_id(name)
                
                valid_fw = [
                    fw for fw, check in fw_checks.items()
                    if check["total"] <= allowed_hits and isoform_id in check["isoforms"]
                ]
                valid_rev = [
                    rv for rv, check in rev_checks.items()
                    if check["total"] <= allowed_hits and isoform_id in check["isoforms"]
                ]
                
                if not valid_fw or not valid_rev:
                    continue
                
                print(f"     Found {len(valid_fw)} FW and {len(valid_rev)} REV primers at this level")
                
                # Try all combinations at this relaxation level
                for fw, rv in product(valid_fw, valid_rev):
                    # Find positions in target sequence
                    fw_pos = seq.find(fw)
                    rv_template = rev_candidates_map[rv]
                    rv_template_pos = seq.find(rv_template)
                    
                    if fw_pos == -1 or rv_template_pos == -1:
                        continue
                    
                    if fw_pos >= rv_template_pos:
                        continue
                    
                    product_size = (rv_template_pos + len(rv_template)) - fw_pos
                    
                    if not (params["prod_min"] <= product_size <= params["prod_max"]):
                        continue
                    
                    # Calculate which isoforms this pair actually amplifies
                    coverage = get_pair_coverage_set(
                        fw, rv, sequences, params["prod_min"], params["prod_max"]
                    )
                    
                    # Ensure target is in coverage (sanity check)
                    if name not in coverage:
                        continue
                    
                    # Accept only if coverage <= allowed_hits
                    if len(coverage) <= allowed_hits:
                        pair_obj = {
                            "forward": fw,
                            "reverse": rv,
                            "product_size": product_size,
                            "coverage": len(coverage),
                            "hits": sorted(list(coverage))
                        }
                        
                        if not found_pair:
                            found_pair = pair_obj
                            print(f"     ✅ Best pair found (amplifies {len(coverage)} isoform(s): {', '.join(coverage)})")
                        
                        if len(found_other_pairs) < 10:
                            found_other_pairs.append(pair_obj)
                
                # If we found primers at this relaxation level, stop trying higher levels
                if found_pair:
                    break
        
        else:  # coverage mode
            # For coverage mode, no relaxation needed - just find valid pairs
            for fw, rv in product(all_fw, all_rev):
                fw_pos = seq.find(fw)
                rv_template = rev_candidates_map[rv]
                rv_template_pos = seq.find(rv_template)
                
                if fw_pos == -1 or rv_template_pos == -1:
                    continue
                
                if fw_pos >= rv_template_pos:
                    continue
                
                product_size = (rv_template_pos + len(rv_template)) - fw_pos
                
                if not (params["prod_min"] <= product_size <= params["prod_max"]):
                    continue
                
                # Coverage mode doesn't need specificity checks
                pair_obj = {
                    "forward": fw,
                    "reverse": rv,
                    "product_size": product_size
                }
                
                if not found_pair:
                    found_pair = pair_obj
                
                if len(found_other_pairs) < 10:
                    found_other_pairs.append(pair_obj)
        
        best_primers[name] = found_pair
        other_options[name] = found_other_pairs if found_other_pairs else []
    
    return {"best": best_primers, "other": other_options}


def get_pair_coverage_set(fw_primer, rv_primer, sequences, prod_min, prod_max):
    """
    Calculate which isoforms a primer pair can amplify.
    
    Args:
        fw_primer: Forward primer sequence (uppercase)
        rv_primer: Reverse primer sequence (uppercase, already RC)
        sequences: List of all sequences
        prod_min: Minimum product size
        prod_max: Maximum product size
    
    Returns:
        Set of isoform IDs that this pair amplifies
    """
    coverage = set()
    
    # rv_primer is already the reverse complement, need template to search
    rv_template = reverse_complement(rv_primer).upper()
    
    for rec in sequences:
        seq = str(rec.seq).upper()
        
        fw_pos = seq.find(fw_primer)
        rv_pos = seq.find(rv_template)
        
        if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
            product_size = (rv_pos + len(rv_template)) - fw_pos
            
            if prod_min <= product_size <= prod_max:
                coverage.add(rec.id)
    
    return coverage

def design_race_primers(sequences, params, mode):
    """
    RACE primer design using windows from sequence ends.
    Specificity mode: Uses progressive relaxation (1 → MAX_RELAX hits) identical to qPCR.
    """
    MAX_RELAX = 4  # Maximum number of isoforms a primer can hit
    
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    
    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()
        
        # Skip sequences too short for window size
        if len(seq) < 2 * params["window_size"]:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Extract terminal regions
        fw_region = seq[:params["window_size"]]
        rev_region = seq[-params["window_size"]:]
        rev_template = reverse_complement(rev_region).upper()
        
        # Find candidate primers in terminal regions
        fw_candidates = find_candidate_primers(
            fw_region, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )
        
        rev_candidates = find_candidate_primers(
            rev_template, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )
        
        if not fw_candidates or not rev_candidates:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Check how many isoforms each primer hits
        fw_checks = {
            fw.upper(): check_cached_primer_hits(fw.upper(), sequences_tuple)
            for fw in fw_candidates
        }
        rev_checks = {
            rv.upper(): check_cached_primer_hits(rv.upper(), sequences_tuple)
            for rv in rev_candidates
        }
        
        # Convert to uppercase
        all_fw = [fw.upper() for fw in fw_candidates]
        all_rev = [rv.upper() for rv in rev_candidates]
        
        found_pair = None
        found_other_pairs = []
        
        if mode == "specificity":
            # PROGRESSIVE RELAXATION for RACE
            for allowed_hits in range(1, MAX_RELAX + 1):
                print(f"  🔍 [{name}] RACE - Trying relaxation level: {allowed_hits} hit(s)...")
                
                isoform_id = extract_isoform_id(name)
                
                valid_fw = [
                    fw for fw, check in fw_checks.items()
                    if check["total"] <= allowed_hits and isoform_id in check["isoforms"]
                ]
                valid_rev = [
                    rv for rv, check in rev_checks.items()
                    if check["total"] <= allowed_hits and isoform_id in check["isoforms"]
                ]
                
                if not valid_fw or not valid_rev:
                    continue
                
                print(f"     Found {len(valid_fw)} FW and {len(valid_rev)} REV primers at this level")
                
                # Try all combinations at this relaxation level
                for fw, rv in product(valid_fw, valid_rev):
                    fw_pos = seq.find(fw)
                    rv_rc = reverse_complement(rv).upper()
                    rv_pos = seq.find(rv_rc)
                    
                    if fw_pos == -1 or rv_pos == -1:
                        continue
                    
                    if fw_pos >= rv_pos:
                        continue
                    
                    product_size = (rv_pos + len(rv_rc)) - fw_pos
                    
                    # Calculate actual coverage for this pair
                    coverage = get_pair_coverage_set_race(
                        fw, rv, sequences
                    )
                    
                    if name not in coverage:
                        continue
                    
                    # Accept only if coverage <= allowed_hits
                    if len(coverage) <= allowed_hits:
                        pair_obj = {
                            "forward": fw,
                            "reverse": rv,
                            "product_size": product_size,
                            "coverage": len(coverage),
                            "hits": sorted(list(coverage))
                        }
                        
                        if not found_pair:
                            found_pair = pair_obj
                            print(f"     ✅ Best pair found (amplifies {len(coverage)} isoform(s): {', '.join(coverage)})")
                        
                        if len(found_other_pairs) < 10:
                            found_other_pairs.append(pair_obj)
                
                if found_pair:
                    break
        
        else:  # coverage mode
            for fw, rv in product(all_fw, all_rev):
                fw_pos = seq.find(fw)
                rv_rc = reverse_complement(rv).upper()
                rv_pos = seq.find(rv_rc)
                
                if fw_pos == -1 or rv_pos == -1:
                    continue
                
                if fw_pos >= rv_pos:
                    continue
                
                product_size = (rv_pos + len(rv_rc)) - fw_pos
                
                pair_obj = {
                    "forward": fw,
                    "reverse": rv,
                    "product_size": product_size
                }
                
                if not found_pair:
                    found_pair = pair_obj
                
                if len(found_other_pairs) < 10:
                    found_other_pairs.append(pair_obj)
        
        best_primers[name] = found_pair
        other_options[name] = found_other_pairs if found_other_pairs else []
    
    return {"best": best_primers, "other": other_options}


def get_pair_coverage_set_race(fw_primer, rv_primer, sequences):
    """
    Calculate which isoforms a RACE primer pair can amplify.
    For RACE, we don't enforce strict product size constraints.
    
    Args:
        fw_primer: Forward primer sequence (uppercase)
        rv_primer: Reverse primer sequence (uppercase, from template strand)
        sequences: List of all sequences
    
    Returns:
        Set of isoform IDs that this pair amplifies
    """
    coverage = set()
    rv_rc = reverse_complement(rv_primer).upper()
    
    for rec in sequences:
        seq = str(rec.seq).upper()
        
        fw_pos = seq.find(fw_primer)
        rv_pos = seq.find(rv_rc)
        
        # For RACE, just check both primers bind in correct orientation
        if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
            coverage.add(rec.id)
    
    return coverage

def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    """
    Finds a minimal set of primer pairs to amplify all isoforms (Set Cover Problem).
    FIXED: Now properly validates qPCR primer pairs with Tm and product size.
    """
    all_isoform_names = {rec.id for rec in sequences}
    params = params_qpcr if primer_type == "qPCR" else params_race
    
    # Pre-compute valid k-mers per isoform with Tm validation
    fw_kmers_by_isoform = {}
    rv_kmers_by_isoform = {}
    
    for rec in sequences:
        seq_str = str(rec.seq)
        
        if primer_type == "qPCR":
            fw_region = seq_str  # Search entire sequence
            rv_region = seq_str
        else:  # RACE
            fw_region = seq_str[:params["window_size"]]
            rv_region = seq_str[-params["window_size"]:]
        
        fw_valid = set()
        rv_valid = set()
        
        # Extract valid forward k-mers with GC AND Tm validation
        for kmer in extract_kmers(fw_region, params["k_range"]):
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if (params["gc_min"] <= gc <= params["gc_max"] and 
                params["tm_min"] <= tm <= params["tm_max"]):
                fw_valid.add(kmer)
        
        # Extract valid reverse k-mers with GC AND Tm validation
        if primer_type == "qPCR":
            rv_region_to_use = rv_region
        else:  # RACE
            rv_region_to_use = reverse_complement(rv_region)
        
        for kmer in extract_kmers(rv_region_to_use, params["k_range"]):
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if (params["gc_min"] <= gc <= params["gc_max"] and 
                params["tm_min"] <= tm <= params["tm_max"]):
                rv_valid.add(kmer)
        
        fw_kmers_by_isoform[rec.id] = fw_valid
        rv_kmers_by_isoform[rec.id] = rv_valid
    
    # Build coverage validation function
    def get_pair_coverage(fw_primer, rv_primer):
        """Returns set of all isoforms that this pair can amplify."""
        covered = set()
        
        for rec in sequences:
            isoform_id = rec.id
            seq_str = str(rec.seq)
            
            # Check if both primers exist for this isoform
            if (fw_primer not in fw_kmers_by_isoform[isoform_id] or 
                rv_primer not in rv_kmers_by_isoform[isoform_id]):
                continue
            
            # Find positions
            fw_pos = seq_str.find(fw_primer)
            
            if primer_type == "qPCR":
                # For qPCR, validate position and product size
                rv_rc = reverse_complement(rv_primer)
                rv_pos = seq_str.find(rv_rc)
                
                if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
                    product_size = (rv_pos + len(rv_rc)) - fw_pos
                    if params["prod_min"] <= product_size <= params["prod_max"]:
                        covered.add(isoform_id)
            else:  # RACE
                # For RACE, just check both primers exist
                if fw_pos != -1:
                    region_3prime = seq_str[-params["window_size"]:]
                    rv_rc = reverse_complement(rv_primer)
                    if rv_rc in region_3prime:
                        covered.add(isoform_id)
        
        return covered
    
    # Generate all possible primer pairs
    all_fw_primers = set()
    all_rv_primers = set()
    for fw_set in fw_kmers_by_isoform.values():
        all_fw_primers.update(fw_set)
    for rv_set in rv_kmers_by_isoform.values():
        all_rv_primers.update(rv_set)
    
    print(f"🔍 Total unique FW primers: {len(all_fw_primers)}")
    print(f"🔍 Total unique RV primers: {len(all_rv_primers)}")
    
    # Create combinations with their coverage
    combinations = []
    for fw in all_fw_primers:
        for rv in all_rv_primers:
            coverage = get_pair_coverage(fw, rv)
            if coverage:
                combinations.append((fw, rv, coverage))
    
    print(f"🔍 Valid primer pair combinations: {len(combinations)}")
    
    # Greedy set cover algorithm
    selected_pairs = []
    pending_isoforms = all_isoform_names.copy()
    
    while pending_isoforms and combinations:
        # Find pair that covers most remaining isoforms
        best_pair = max(combinations, key=lambda item: len(item[2].intersection(pending_isoforms)))
        fw, rv, covered_by_best = best_pair
        
        newly_covered = covered_by_best.intersection(pending_isoforms)
        
        if not newly_covered:
            break
        
        # Store ALL isoforms this pair amplifies (not just newly covered)
        selected_pairs.append({
            "forward": fw,
            "reverse": rv,
            "isoforms": sorted(list(covered_by_best))
        })
        
        print(f"✅ Selected pair covers {len(newly_covered)} new isoforms (total: {len(covered_by_best)})")
        
        pending_isoforms -= newly_covered
        combinations.remove(best_pair)
    
    if pending_isoforms:
        print(f"⚠️ Warning: Could not cover {len(pending_isoforms)} isoforms: {pending_isoforms}")
    
    return {"selected": selected_pairs, "pending": list(pending_isoforms)}

# ============================================================================== 
# --- SCRIPT EXECUTION --- 
# ==============================================================================
def parse_params(qpcr_params, race_params):
    if qpcr_params:
        qpcr_params = {
            "k_range": range(
                qpcr_params["k_range"][0], 
                qpcr_params["k_range"][1]
            ),   
            "gc_min": qpcr_params["gc_range"][0]*100,
            "gc_max": qpcr_params["gc_range"][1]*100,
            "prod_min": qpcr_params["prod_range"][0],
            "prod_max": qpcr_params["prod_range"][1],
            "tm_min": qpcr_params["tm_range"][0],
            "tm_max": qpcr_params["tm_range"][1],
        }
    
    if race_params:
        race_params = {
            "window_size": race_params["window_size"],
            "k_range": range(
                race_params["k_range"][0], 
                race_params["k_range"][1]
            ),   
            "gc_min": race_params["gc_range"][0]*100,
            "gc_max": race_params["gc_range"][1]*100,
            "tm_min": race_params["tm_range"][0],
            "tm_max": race_params["tm_range"][1],
        }
    
    return qpcr_params, race_params

def create_design(fasta_file_path, primer_type, design_mode, qpcr_params, race_params):
    """
    Main function to execute the primer design pipeline.
    Returns:
        - For specificity mode: tuple of (num_sequences, best_pairs_df, other_pairs_df)
        - For coverage mode: tuple of (num_sequences, coverage_df, None)
    """
    # Clear cache at start of each run for accurate statistics
    clear_tm_cache()
    
    print(f"📂 Reading FASTA file: {fasta_file_path}")
    num_sequences = 0
    try:
        sequence_set = list(SeqIO.parse(fasta_file_path, "fasta"))
        if not sequence_set:
            print("❌ Error: No sequences found in the FASTA file.")
            return None, None, None
        num_sequences = len(sequence_set)
        print(f"✅ Loaded {num_sequences} sequences.")
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at '{fasta_file_path}'")
        return None, None, None
    
    qpcr_params, race_params = parse_params(qpcr_params, race_params)
    results = None
    
    if design_mode == "specificity":
        print(f"🧬 Running in SPECIFICITY mode for {primer_type} primers...")
        if primer_type == "qPCR":
            results = design_qpcr_primers(sequence_set, qpcr_params, design_mode)
        elif primer_type == "RACE":
            results = design_race_primers(sequence_set, race_params, design_mode)
        else:
            print(f"⚠️ Unsupported primer type: {primer_type}")
            return None, None, None
    
    elif design_mode == "coverage":
        print(f"🎯 Running in COVERAGE mode for {primer_type} primers...")
        results = find_coverage_primers(sequence_set, qpcr_params, race_params, primer_type)
    
    else:
        print(f"⚠️ Unsupported design mode: {design_mode}")
        return None, None, None
    
    # Print cache statistics for performance monitoring
    cache_stats = get_cache_stats()
    print(f"\n⚡ Performance Stats:")
    print(f"   Tm Cache: {cache_stats['cache_size']} unique sequences")
    print(f"   Cache hits: {cache_stats['hits']} ({cache_stats['hit_rate']:.1f}%)")
    print(f"   Cache saved ~{cache_stats['hits']} Tm calculations!")
    
    # --- Process and Return Results ---
    if design_mode == "specificity":
        best_df_rows = []
        other_df_rows = []
        
        for isoform_name, pair in results["best"].items():
            if not pair: 
                continue
            f, r = pair["forward"], pair["reverse"]
            best_df_rows.append({
                "Isoform": isoform_name,
                "Primer_Forward": f,
                "Primer_Reverse": r,
                "GC_Forward": calculate_gc_content(f),
                "GC_Reverse": calculate_gc_content(r),
                "Tm_Forward": calculate_melting_temp(f),
                "Tm_Reverse": calculate_melting_temp(r),
                "Hairpin_Forward": find_longest_complementary_run(f, reverse_complement(f)),
                "Hairpin_Reverse": find_longest_complementary_run(r, reverse_complement(r)),
                "SelfDimer_Forward": find_longest_complementary_run(f, f),
                "SelfDimer_Reverse": find_longest_complementary_run(r, r),
                "CrossDimer": find_longest_complementary_run(f, reverse_complement(r)),
                "Product_Size": pair.get("product_size", None)
            })
        
        for isoform_name, pairs in results["other"].items():
            for pair in pairs:
                if not pair:
                    continue
                f, r = pair["forward"], pair["reverse"]
                other_df_rows.append({
                    "Isoform": isoform_name,
                    "Primer_Forward": f,
                    "Primer_Reverse": r,
                    "GC_Forward": calculate_gc_content(f),
                    "GC_Reverse": calculate_gc_content(r),
                    "Tm_Forward": calculate_melting_temp(f),
                    "Tm_Reverse": calculate_melting_temp(r),
                    "Hairpin_Forward": find_longest_complementary_run(f, reverse_complement(f)),
                    "Hairpin_Reverse": find_longest_complementary_run(r, reverse_complement(r)),
                    "SelfDimer_Forward": find_longest_complementary_run(f, f),
                    "SelfDimer_Reverse": find_longest_complementary_run(r, r),
                    "CrossDimer": find_longest_complementary_run(f, reverse_complement(r)),
                    "Product_Size": pair.get("product_size", None)
                })
        
        best_df = pd.DataFrame(best_df_rows) if best_df_rows else None
        other_df = pd.DataFrame(other_df_rows) if other_df_rows else None
        
        if best_df is not None and not best_df.empty: 
            best_df.index = range(1, len(best_df) + 1)
        if other_df is not None and not other_df.empty: 
            other_df.index = range(1, len(other_df) + 1)
        
        print("✅ Finished specificity design.")
        return num_sequences, best_df, other_df
    
    elif design_mode == "coverage":
        coverage_df_rows = []
        for idx, item in enumerate(results["selected"], 1):
            f, r = item["forward"], item["reverse"]
            coverage_df_rows.append({
                "Pair_ID": f"Pair_{idx}",
                "Primer_Forward": f,
                "Primer_Reverse": r,
                "GC_Forward": calculate_gc_content(f),
                "GC_Reverse": calculate_gc_content(r),
                "Tm_Forward": calculate_melting_temp(f),
                "Tm_Reverse": calculate_melting_temp(r),
                "Num_Isoforms": len(item["isoforms"]),
                "Covered_Isoforms": ";".join(item["isoforms"]),
            })
        
        coverage_df = pd.DataFrame(coverage_df_rows) if coverage_df_rows else None
        if coverage_df is not None and not coverage_df.empty: 
            coverage_df.index = range(1, len(coverage_df) + 1)
        
        print(f"✅ Finished coverage design.")
        print(f"📊 Found {len(results['selected'])} primer pairs")
        print(f"📊 Uncovered isoforms: {len(results['pending'])}")
        if results['pending']:
            print(f"⚠️ Uncovered: {', '.join(results['pending'])}")
        
        return num_sequences, coverage_df, None

if __name__ == "__main__":
    # testing (locally)
    FASTA_FILE_PATH = "./sample_data/avengers-3.fa"
    PRIMER_TYPE = "RACE"
    DESIGN_MODE = "specificity" 
    PARAMS_QPCR = {
        'k_range': (18, 26), 
        'gc_range': (0.5, 0.7), 
        'prod_range': (80, 150), 
        'tm_range': (57.0, 62.0)
    }
    PARAMS_RACE = {
        'window_size': 200, 
        'k_range': (23, 29), 
        'gc_range': (0.5, 0.7), 
        'tm_range': (57.0, 62.0)
    }
    create_design(FASTA_FILE_PATH, PRIMER_TYPE, DESIGN_MODE, PARAMS_QPCR, PARAMS_RACE)
