# ============================================================================== 
# primer_designer.py
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
    Calculates the GC content of a DNA sequence (percentage 0-100).
    Works regardless of case by normalizing to upper.
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc_count = seq.count("G") + seq.count("C")
    return round(gc_count / len(seq) * 100, 2)


def extract_kmers(sequence, k_range):
    """
    Yields all k-mers of specified lengths from a sequence (uppercased).
    """
    if not sequence:
        return
    seq = sequence.upper()
    for k in k_range:
        if 0 < k <= len(seq):
            for i in range(len(seq) - k + 1):
                yield seq[i:i+k]


COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def are_bases_complementary(base1, base2):
    """Checks if two DNA bases are complementary using a fast lookup."""
    return COMPLEMENT_MAP.get(base1) == base2

def reverse_complement(sequence):
    """Computes the reverse complement of a DNA sequence using Biopython (returns UPPER)."""
    return str(Seq(sequence.upper()).reverse_complement())

def find_longest_complementary_run(seq1, seq2):
    """
    Finds the length of the longest perfectly complementary substring between two sequences.
    Memory-optimized DP approach.
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

def calculate_melting_temp(sequence, na_conc=50e-3):
    """
    Calculates the melting temperature (Tm) using Biopython's nearest-neighbor model.
    Normalizes input to string/upper.
    """
    try:
        return round(mt.Tm_NN(str(sequence).upper(), Na=na_conc), 1)
    except Exception:
        return 0.0

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
    Finds all sequence records that contain a given primer sequence or its reverse complement.
    Returns a dict: {"total": int, "isoforms": [id1, id2, ...]}.
    This function normalizes everything to uppercase for exact matching.
    """
    p = primer.upper()
    primer_rc = reverse_complement(p).upper()
    hits = [
        rec.id
        for rec in sequences
        if p in str(rec.seq).upper() or primer_rc in str(rec.seq).upper()
    ]
    return {"total": len(hits), "isoforms": hits}


def find_candidate_primers(sequence, k_range, gc_min, gc_max, tm_min, tm_max):
    """
    Return candidate k-mers (UPPERCASE) that meet GC and Tm constraints.
    Uses calculate_melting_temp for consistency.
    """
    candidates = []
    seq = sequence.upper()
    for k in k_range:
        if k <= 0 or k > len(seq):
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
                candidates.append(kmer)
    return candidates

# ============================================================================== 
# --- PRIMER DESIGN CORE LOGIC --- 
# ==============================================================================

@lru_cache(maxsize=10000)
def check_cached_primer_hits(primer, sequences_tuple):
    """
    Cached wrapper for check_primer_hits.
    NOTE: callers should pass primer in UPPERCASE to maximize cache hits.
    sequences_tuple should contain (id, seq) with seq uppercased (prepare_sequences_for_caching does that).
    """
    # sequences_tuple already produced by prepare_sequences_for_caching (uppercased seqs)
    sequences = [SeqIO.SeqRecord(Seq(seq), id=name) for name, seq in sequences_tuple]
    return check_primer_hits(primer, sequences)


def prepare_sequences_for_caching(sequences):
    """
    Converts a list of SeqRecord objects to a hashable tuple for caching.
    Also normalizes sequences to UPPERCASE here.
    """
    return tuple((rec.id, str(rec.seq).upper()) for rec in sequences)


def check_pair_off_target(fw_primer, rv_primer, target_isoform, sequences, prod_min, prod_max):
    """
    Check if a primer pair amplifies any isoform other than the target.
    fw_primer: forward primer (UPPER)
    rv_primer: actual reverse primer sequence (UPPER, i.e., reverse complement)
    target_isoform: ID of the target isoform (string)
    sequences: list of SeqRecord
    """
    # compute the template for reverse search (reverse of rv_primer)
    rv_template = reverse_complement(rv_primer).upper()
    fw = fw_primer.upper()
    for rec in sequences:
        if rec.id == target_isoform:
            continue
        other_seq = str(rec.seq).upper()
        fw_pos = other_seq.find(fw)
        rv_pos = other_seq.find(rv_template)
        if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
            other_product_size = (rv_pos + len(rv_template)) - fw_pos
            if prod_min <= other_product_size <= prod_max:
                # Off-target found
                # debug print
                print(f"  ⚠️ Off-target detected: {rec.id} (product size: {other_product_size} bp)")
                return True
    return False


def design_qpcr_primers(sequences, params, mode):
    """
    Main logic for designing qPCR primers for either specificity or coverage.
    For specificity mode: ensures primers are unique AND pairs don't amplify off-target isoforms.
    """
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
        
        # Candidate k-mers from full sequence (returned uppercase)
        fw_candidates = find_candidate_primers(
            seq,
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        rev_template_candidates = find_candidate_primers(
            seq,
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        
        # Map reverse complement (actual primer) -> template k-mer (all uppercase)
        rev_candidates_map = {reverse_complement(p).upper(): p.upper() for p in rev_template_candidates}
        
        if not fw_candidates or not rev_candidates_map:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Which isoforms each primer hits (use cached checks; pass UPPER primers)
        fw_checks = {fw.upper(): check_cached_primer_hits(fw.upper(), sequences_tuple) for fw in fw_candidates}
        rev_checks = {rv.upper(): check_cached_primer_hits(rv.upper(), sequences_tuple) for rv in rev_candidates_map.keys()}
        
        valid_fw = [fw.upper() for fw in fw_candidates]
        valid_rev = list(rev_candidates_map.keys())  # already upper
        
        if mode == "specificity":
            isoform_id = extract_isoform_id(name)
            
            # require primer to be unique (only hits this isoform)
            valid_fw = [
                fw for fw, chk in fw_checks.items()
                if chk["total"] == 1 and isoform_id in chk["isoforms"]
            ]
            valid_rev = [
                rv for rv, chk in rev_checks.items()
                if chk["total"] == 1 and isoform_id in chk["isoforms"]
            ]
        
        print(f"🔍 [{name}] Found {len(valid_fw)} unique FW primers, {len(valid_rev)} unique REV primers.")
        
        best_pair = None
        other_pairs = []
        
        for fw, rv in product(valid_fw, valid_rev):
            # fw, rv are uppercase primers (rv is actual primer seq - reverse complement)
            fw_pos = seq.find(fw)
            rv_template = rev_candidates_map[rv]  # template (uppercase)
            rv_template_pos = seq.find(rv_template)
            
            if fw_pos == -1 or rv_template_pos == -1:
                continue
            
            # Ensure forward is before reverse
            if fw_pos >= rv_template_pos:
                continue
            
            # Calculate product size (template-based)
            product_size = (rv_template_pos + len(rv_template)) - fw_pos
            
            if not (params["prod_min"] <= product_size <= params["prod_max"]):
                continue
            
            # Pair-level specificity check: ensure no other isoform produces an amplicon
            if mode == "specificity":
                off_target = check_pair_off_target(
                    fw, rv, name, sequences, params["prod_min"], params["prod_max"]
                )
                if off_target:
                    continue  # skip pair that amplifies another isoform
            
            primer_pair = {
                "forward": fw,
                "reverse": rv,
                "product_size": product_size
            }
            
            if mode == "specificity":
                if not best_pair:
                    best_pair = primer_pair
                if len(other_pairs) < 10:
                    other_pairs.append(primer_pair)
        
        best_primers[name] = best_pair
        other_options[name] = other_pairs
    
    return {"best": best_primers, "other": other_options}


def design_race_primers(sequences, params, mode):
    """
    Main logic for designing RACE primers for either specificity or coverage.
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    
    valid_sequences = [(rec.id, str(rec.seq)) for rec in sequences
                       if len(str(rec.seq)) >= 2 * params["window_size"]]
    
    if not valid_sequences:
        return {"best": best_primers, "other": other_options}
    
    for name, seq in valid_sequences:
        seq = seq.upper()
        fw_region = seq[:params["window_size"]]
        rev_region = seq[-params["window_size"]:]
        rev_template = reverse_complement(rev_region)  # uppercase
        
        fw_candidates = find_candidate_primers(
            fw_region,
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        
        rev_candidates = find_candidate_primers(
            rev_template,
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        
        print(f"\n🧬 Isoform: {name}")
        print(f"Forward region: {fw_region[:50]}...")
        print(f"Reverse region: {rev_region[:50]}...")
        print(f"Rev template strand: {rev_template[:50]}...")
        print(f"First FW primers: {fw_candidates[:3]}")
        print(f"First REV primers: {rev_candidates[:3]}")
        
        if not fw_candidates or not rev_candidates:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        fw_checks = {fw: check_cached_primer_hits(fw, sequences_tuple) for fw in fw_candidates}
        rev_checks = {rv: check_cached_primer_hits(rv, sequences_tuple) for rv in rev_candidates}
        
        valid_fw = fw_candidates
        valid_rev = rev_candidates
        
        if mode == "specificity":
            isoform_id = extract_isoform_id(name)
            valid_fw = [fw for fw, check in fw_checks.items() if isoform_id in check["isoforms"]]
            valid_rev = [rv for rv, check in rev_checks.items() if isoform_id in check["isoforms"]]
        
        print(f"Valid FW (specific): {[p[:10] for p in valid_fw]}")
        print(f"Valid REV (specific): {[p[:10] for p in valid_rev]}")
        
        best_pair = None
        other_pairs = []
        
        for fw, rv in product(valid_fw, valid_rev):
            primer_pair = {"forward": fw, "reverse": rv, "product_size": None}
            if mode == "specificity":
                if not best_pair:
                    best_pair = primer_pair
                if len(other_pairs) < 10:
                    other_pairs.append(primer_pair)
        
        best_primers[name] = best_pair
        other_options[name] = other_pairs
    
    return {"best": best_primers, "other": other_options}


def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    """
    Finds a minimal set of primer pairs to amplify all isoforms (Set Cover Problem).
    The implementation normalizes sequences to uppercase and validates Tm/product size.
    Note: for very large datasets this may be slow (combinatorial).
    """
    all_isoform_names = {rec.id for rec in sequences}
    params = params_qpcr if primer_type == "qPCR" else params_race
    
    # Pre-compute valid k-mers per isoform with Tm validation (uppercase)
    fw_kmers_by_isoform = {}
    rv_kmers_by_isoform = {}
    
    for rec in sequences:
        seq_str = str(rec.seq).upper()
        
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
    
    # Build all valid primer pairs and their coverage (may be large)
    def get_pair_coverage(fw_primer, rv_primer):
        """Returns set of all isoforms that this pair can amplify."""
        covered = set()
        for rec in sequences:
            isoform_id = rec.id
            seq_str = str(rec.seq).upper()
            
            # Check if both primers exist for this isoform
            if (fw_primer not in fw_kmers_by_isoform[isoform_id] or 
                rv_primer not in rv_kmers_by_isoform[isoform_id]):
                continue
            
            # Find positions
            fw_pos = seq_str.find(fw_primer)
            
            if primer_type == "qPCR":
                # For qPCR, reverse primer is the actual primer (not template)
                rv_rc = reverse_complement(rv_primer)
                rv_pos = seq_str.find(rv_rc)
                
                if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
                    product_size = (rv_pos + len(rv_rc)) - fw_pos
                    if params["prod_min"] <= product_size <= params["prod_max"]:
                        covered.add(isoform_id)
            else:  # RACE
                # For RACE, just check both primers exist in their regions
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
# --- SCRIPT EXECUTION / PARAM PARSING --- 
# ==============================================================================

def parse_params(qpcr_params, race_params):
    """
    Parse UI param shapes (qpcr: k_range,gc_range,prod_range,tm_range).
    Returns normalized dicts with ranges and numeric thresholds.
    """
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
    # quick local test (adjust FASTA path if needed)
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
