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
def calculate_melting_temp(sequence, na_conc=50e-3):
    """
    Calculates the melting temperature (Tm) using Biopython's nearest-neighbor model.
    """
    try:
        return round(mt.Tm_NN(str(sequence), Na=na_conc), 1)
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
    candidates = []
    for k in k_range:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            gc = calculate_gc_content(kmer)
            tm = mt.Tm_NN(Seq(kmer))
            if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
                candidates.append(kmer)
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
    Main logic for designing qPCR primers for either specificity or coverage.
    NOTE: This version scans the ENTIRE sequence for qPCR primers (requested change).
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    for rec in sequences:
        name = rec.id
        seq = str(rec.seq)
        if len(seq) < params["prod_min"]:
            best_primers[name] = None
            other_options[name] = []
            continue
        # Scan the full sequence for forward and reverse-template candidates
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
        # Map reverse complement (actual primer) -> template k-mer
        rev_candidates_map = {reverse_complement(p): p for p in rev_template_candidates}
        if not fw_candidates or not rev_candidates_map:
            best_primers[name] = None
            other_options[name] = []
            continue
        fw_checks = {fw: check_cached_primer_hits(fw, sequences_tuple) for fw in fw_candidates}
        rev_checks = {rv: check_cached_primer_hits(rv, sequences_tuple) for rv in rev_candidates_map.keys()}
        valid_fw = fw_candidates
        valid_rev = list(rev_candidates_map.keys())
        if mode == "specificity":
            isoform_id = extract_isoform_id(name)
            valid_fw = [
                fw for fw, check in fw_checks.items()
                if isoform_id in check["isoforms"]
            ]
            valid_rev = [
                rv for rv, check in rev_checks.items()
                if isoform_id in check["isoforms"]
            ]
        print(f"🔍 [{name}] Found {len(valid_fw)} FW primers, {len(valid_rev)} REV primers.")
        best_pair = None
        other_pairs = []
        for fw, rv in product(valid_fw, valid_rev):
            # Find positions in sequence
            fw_pos = seq.find(fw)
            rv_template = rev_candidates_map[rv]
            rv_template_pos = seq.find(rv_template)
            if fw_pos == -1 or rv_template_pos == -1:
                continue
            # Ensure forward is before reverse
            if fw_pos >= rv_template_pos:
                continue
            # Calculate product size
            product_size = (rv_template_pos + len(rv_template)) - fw_pos
            if not (params["prod_min"] <= product_size <= params["prod_max"]):
                continue
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
        fw_region = seq[:params["window_size"]]
        rev_region = seq[-params["window_size"]:]
        rev_template = reverse_complement(rev_region)
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
    Greedy set-cover for primer pairs.

    For qPCR: scans full sequences but avoids enumerating all fw×rv combinations.
    For RACE: uses windows at 5' and 3' ends (existing behavior).
    """
    # Map id -> sequence string (uppercase)
    seq_map = {rec.id: str(rec.seq).upper() for rec in sequences}
    all_isoform_names = set(seq_map.keys())
    params = params_qpcr if primer_type == "qPCR" else params_race

    # Helper to get regions
    def get_regions(seq_str):
        if primer_type == "qPCR":
            return seq_str, seq_str
        else:
            return seq_str[: params["window_size"]], seq_str[-params["window_size"] :]

    # Build fw_kmers_by_isoform and rv_kmers_by_isoform (store kmers UPPER)
    fw_kmers_by_isoform = {}
    rv_kmers_by_isoform = {}
    for iso, seq_str in seq_map.items():
        fw_region, rv_region = get_regions(seq_str)
        fw_set = set(k.upper() for k in extract_kmers(fw_region, params["k_range"]))
        rv_set = set(k.upper() for k in extract_kmers(rv_region, params["k_range"]))
        # Filter by GC and Tm
        fw_valid = {k for k in fw_set if params["gc_min"] <= calculate_gc_content(k) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(k) <= params["tm_max"]}
        rv_valid = {k for k in rv_set if params["gc_min"] <= calculate_gc_content(k) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(k) <= params["tm_max"]}
        fw_kmers_by_isoform[iso] = fw_valid
        rv_kmers_by_isoform[iso] = rv_valid

    # Build inverted indexes: primer -> isoforms it appears in
    fw_index = {}
    rv_index = {}
    for iso in seq_map:
        for k in fw_kmers_by_isoform.get(iso, set()):
            fw_index.setdefault(k, set()).add(iso)
        for k in rv_kmers_by_isoform.get(iso, set()):
            rv_index.setdefault(k, set()).add(iso)

    # Prepare candidate lists
    all_fw = list(fw_index.keys())
    all_rv = list(rv_index.keys())

    # Greedy selection
    selected_pairs = []
    pending = set(all_isoform_names)

    # Quick bailout: no candidates
    if not all_fw or not all_rv:
        return {"selected": [], "pending": sorted(list(pending))}

    # For qPCR we must check product size constraints; for RACE we do simpler coverage check
    if primer_type == "qPCR":
        prod_min = params["prod_min"]
        prod_max = params["prod_max"]

        # Precompute for speed reverse_complements of rv templates
        rv_to_rc = {rv: reverse_complement(rv) for rv in all_rv}

        # Greedy loop: instead of enumerating all fw×rv, evaluate promising fw and for each find best rv
        while pending:
            best_pair = None
            best_cover = set()

            # Optionally limit number of fw candidates considered per round if too many
            # e.g., consider fw candidates sorted by how many pending isoforms they hit
            fw_candidates_sorted = sorted(all_fw, key=lambda f: len(fw_index.get(f, set()).intersection(pending)), reverse=True)

            # limit to top N forwards to speed up (tune as needed). Keep large for accuracy.
            MAX_FW_TO_TEST = min(len(fw_candidates_sorted), 200)
            for fw in fw_candidates_sorted[:MAX_FW_TO_TEST]:
                # isoforms where this fw exists and are still pending
                iso_candidates = fw_index.get(fw, set()).intersection(pending)
                if not iso_candidates:
                    continue

                # For this forward, build a map rv -> set of isoforms (subset of iso_candidates) it can pair with
                rv_cov_map = {}

                for iso in iso_candidates:
                    seq = seq_map[iso]
                    fw_pos = seq.find(fw)
                    if fw_pos == -1:
                        continue
                    # For each reverse-template present in this isoform, check if its RC pairs with fw (positions, product size)
                    for rv_template in rv_kmers_by_isoform.get(iso, set()):
                        rv_rc = rv_to_rc.get(rv_template)
                        if not rv_rc:
                            # if rv_template not in rv_to_rc keys (shouldn't happen), compute:
                            rv_rc = reverse_complement(rv_template)
                            rv_to_rc[rv_template] = rv_rc
                        rv_pos = seq.find(rv_rc)
                        if rv_pos == -1:
                            continue
                        if fw_pos >= rv_pos:  # forward must be upstream of reverse
                            continue
                        product_size = (rv_pos + len(rv_rc)) - fw_pos
                        if not (prod_min <= product_size <= prod_max):
                            continue
                        rv_cov_map.setdefault(rv_rc, set()).add(iso)

                # choose best rv for this fw
                if not rv_cov_map:
                    continue
                # pick rv that covers most isoforms (intersection with pending)
                best_rv_for_fw, covered_iso = max(rv_cov_map.items(), key=lambda it: len(it[1].intersection(pending)))
                covered_iso = covered_iso.intersection(pending)
                if len(covered_iso) > len(best_cover):
                    best_cover = covered_iso
                    best_pair = (fw, best_rv_for_fw)

            # If no pair adds coverage, break
            if not best_pair or not best_cover:
                break

            # record pair and remove covered isoforms
            selected_pairs.append({
                "forward": best_pair[0],
                "reverse": best_pair[1],
                "isoforms": sorted(list(best_cover))
            })
            pending -= best_cover

        return {"selected": selected_pairs, "pending": sorted(list(pending))}

    else:
        # RACE: simpler greedy using intersections of k-mer hit sets (unchanged)
        combinations = []
        for fw_kmer, fw_hits in fw_index.items():
            for rv_kmer, rv_hits in rv_index.items():
                common_hits = fw_hits.intersection(rv_hits)
                if common_hits:
                    combinations.append((fw_kmer, rv_kmer, common_hits))

        selected = []
        while pending and combinations:
            best_pair = max(combinations, key=lambda item: len(item[2].intersection(pending)))
            fw, rv, covered_by_best = best_pair
            newly_covered = covered_by_best.intersection(pending)
            if not newly_covered:
                break
            selected.append({"forward": fw, "reverse": rv, "isoforms": list(newly_covered)})
            pending -= newly_covered
            combinations.remove(best_pair)

        return {"selected": selected, "pending": sorted(list(pending))}
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
    print(f"📂 Reading FASTA file: {fasta_file_path}")
    num_sequences = 0
    try:
        sequence_set = list(SeqIO.parse(fasta_file_path, "fasta"))
        if not sequence_set:
            print("❌ Error: No sequences found in the FASTA file.")
            return (None, None) if design_mode == "specificity" else None
        num_sequences = len(sequence_set)
        print(f"✅ Loaded {num_sequences} sequences.")
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at '{fasta_file_path}'")
        return (None, None) if design_mode == "specificity" else None
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
            return None, None
    elif design_mode == "coverage":
        print(f"🎯 Running in COVERAGE mode for {primer_type} primers...")
        results = find_coverage_primers(sequence_set, qpcr_params, race_params, primer_type)
    else:
        print(f"⚠️ Unsupported design mode: {design_mode}")
        return (None, None) if design_mode == "specificity" else None
    # --- Process and Return Results ---
    if design_mode == "specificity":
        # Process best pairs
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
        # Process other pairs (up to 10 per isoform)
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
        if best_df is not None and not best_df.empty: best_df.index = range(1, len(best_df) + 1)
        if other_df is not None and not other_df.empty: other_df.index = range(1, len(other_df) + 1)
        
        print("✅ Finished specificity design.")
        return num_sequences, best_df, other_df
    elif design_mode == "coverage":
        coverage_df_rows = []
        for item in results["selected"]:
            coverage_df_rows.append({
                "Primer_Forward": item["forward"],
                "Primer_Reverse": item["reverse"],
                "Covered_Isoforms": ";".join(item["isoforms"]),
            })
        coverage_df = pd.DataFrame(coverage_df_rows) if coverage_df_rows else None
        if coverage_df is not None and not coverage_df.empty: coverage_df.index = range(1, len(coverage_df) + 1)
        print("✅ Finished coverage design.")
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
