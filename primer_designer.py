# ============================================================================== 
# --- IMPORTS --- 
# ============================================================================== 
import re
from itertools import product
from functools import lru_cache
from bisect import bisect_left, bisect_right
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# ============================================================================== 
# --- CORE BIOINFORMATICS UTILITY FUNCTIONS --- 
# ============================================================================== 
def calculate_gc_content(sequence):
    if not sequence:
        return 0.0
    seq = str(sequence).upper()
    gc_count = seq.count("G") + seq.count("C")
    return round(gc_count / len(seq) * 100, 2)

def extract_kmers(sequence, k_range):
    if not sequence:
        return
    seq = str(sequence).upper()
    for k in k_range:
        if 0 < k <= len(seq):
            for i in range(len(seq) - k + 1):
                yield seq[i:i+k]

COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}

def are_bases_complementary(base1, base2):
    return COMPLEMENT_MAP.get(base1) == base2

def reverse_complement(sequence):
    return str(Seq(str(sequence)).reverse_complement())

def find_longest_complementary_run(seq1, seq2):
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
    try:
        return round(mt.Tm_NN(str(sequence), Na=na_conc), 1)
    except Exception:
        return 0.0

def extract_isoform_id(name):
    match = re.search(r"CASC15-(\d+)", name)
    if not match:
        print(f"⚠️ Could not extract isoform ID from: {name}")
    return match.group(1) if match else ""

def check_primer_hits(primer, sequences):
    primer_uc = str(primer).upper()
    primer_rc = reverse_complement(primer_uc)
    hits = [
        rec.id
        for rec in sequences
        if primer_uc in str(rec.seq).upper() or primer_rc in str(rec.seq).upper()
    ]
    return {"total": len(hits), "isoforms": ",".join(hits)}

# ============================================================================== 
# --- PRIMER FINDERS --- 
# ============================================================================== 
def find_candidate_primers(sequence, k_range, gc_min, gc_max, tm_min, tm_max):
    seq = str(sequence).upper()
    candidates = []
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

def find_candidate_primers_with_positions(sequence, k_range, gc_min, gc_max, tm_min, tm_max):
    seq = str(sequence).upper()
    candidates = []
    for k in k_range:
        if k <= 0 or k > len(seq):
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
                candidates.append((kmer, i))
    return candidates

# ============================================================================== 
# --- PRIMER DESIGN CORE LOGIC --- 
# ============================================================================== 
@lru_cache(maxsize=10000)
def check_cached_primer_hits(primer, sequences_tuple):
    sequences = [SeqIO.SeqRecord(Seq(seq), id=name) for name, seq in sequences_tuple]
    return check_primer_hits(primer, sequences)

def prepare_sequences_for_caching(sequences):
    return tuple((rec.id, str(rec.seq)) for rec in sequences)

# ============================================================================== 
# PAIR COVERAGE FUNCTIONS
# ============================================================================== 
def get_pair_coverage_set(fw_primer, rv_primer, sequences, prod_min, prod_max):
    coverage = set()
    rv_template = reverse_complement(rv_primer)
    for rec in sequences:
        seq = str(rec.seq).upper()
        fw_pos = seq.find(fw_primer.upper())
        rv_pos = seq.find(rv_template.upper())
        if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
            product_size = (rv_pos + len(rv_template)) - fw_pos
            if prod_min <= product_size <= prod_max:
                coverage.add(rec.id)
    return coverage

def get_pair_coverage_set_race(fw_primer, rv_primer, sequences):
    coverage = set()
    for rec in sequences:
        seq = str(rec.seq).upper()
        fw_pos = seq.find(fw_primer.upper())
        rv_pos = seq.find(rv_primer.upper())
        if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
            coverage.add(rec.id)
    return coverage

# ============================================================================== 
# MAIN DESIGN FUNCTIONS - qPCR (forward-first, position-based) 
# ============================================================================== 
def design_qpcr_primers(sequences, params, mode):
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    MAX_RELAX = 4
    tm_mid = (params["tm_min"] + params["tm_max"]) / 2.0

    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()

        if len(seq) < params["prod_min"]:
            best_primers[name] = None
            other_options[name] = []
            continue

        fw_pos_cands_all = find_candidate_primers_with_positions(
            seq, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )
        rev_pos_cands_all = find_candidate_primers_with_positions(
            seq, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )

        if not fw_pos_cands_all or not rev_pos_cands_all:
            best_primers[name] = None
            other_options[name] = []
            continue

        rev_list_by_pos = sorted(rev_pos_cands_all, key=lambda kv: kv[1])
        rev_positions = [pos for (_kmer, pos) in rev_list_by_pos]

        pair_candidates = []

        for fw_kmer, fw_pos in fw_pos_cands_all:
            rv_search_min = fw_pos + 1
            rv_search_max = fw_pos + params["prod_max"]

            left_idx = bisect_left(rev_positions, rv_search_min)
            right_idx = bisect_right(rev_positions, rv_search_max)

            for rv_kmer, rv_pos in rev_list_by_pos[left_idx:right_idx]:
                if fw_pos >= rv_pos:
                    continue

                L_rv = len(rv_kmer)
                product_size = (rv_pos + L_rv) - fw_pos
                if not (params["prod_min"] <= product_size <= params["prod_max"]):
                    continue

                rv_actual = reverse_complement(rv_kmer)

                coverage = get_pair_coverage_set(fw_kmer, rv_actual, sequences, params["prod_min"], params["prod_max"])
                if name not in coverage:
                    continue

                off_targets = len(coverage) - 1
                tm_fw = calculate_melting_temp(fw_kmer)
                tm_rv = calculate_melting_temp(rv_actual)
                tm_diff = abs(tm_fw - tm_rv)
                cross_dimer = find_longest_complementary_run(fw_kmer, reverse_complement(rv_actual))
                prod_mid = (params["prod_min"] + params["prod_max"]) / 2.0
                prod_dist = abs(product_size - prod_mid)

                pair_candidates.append({
                    "forward": fw_kmer,
                    "fw_pos": fw_pos,
                    "reverse": rv_actual,
                    "rv_pos": rv_pos,
                    "product_size": product_size,
                    "coverage": coverage,
                    "off_targets": off_targets,
                    "tm_diff": tm_diff,
                    "cross_dimer": cross_dimer,
                    "prod_dist": prod_dist
                })

        print(f"🔍 [{name}] Generated {len(pair_candidates)} position-based pair candidates")

        if not pair_candidates:
            best_primers[name] = None
            other_options[name] = []
            continue

        if mode == "specificity":
            found_pair = None
            other_pairs = []
            for allowed_off in range(0, MAX_RELAX + 1):
                filtered = [p for p in pair_candidates if p["off_targets"] <= allowed_off]
                if not filtered:
                    continue
                print(f"     🎯 Relaxation level {allowed_off}: {len(filtered)} pairs available")
                filtered_sorted = sorted(
                    filtered,
                    key=lambda p: (p["off_targets"], p["tm_diff"], p["cross_dimer"], p["prod_dist"])
                )
                best = filtered_sorted[0]
                found_pair = {
                    "forward": best["forward"],
                    "reverse": best["reverse"],
                    "product_size": best["product_size"]
                }
                for alt in filtered_sorted[:10]:
                    other_pairs.append({
                        "forward": alt["forward"],
                        "reverse": alt["reverse"],
                        "product_size": alt["product_size"]
                    })
                print(f"     ✅ Best pair: {best['off_targets']} off-targets, Tm diff={best['tm_diff']:.1f}°C")
                break
            best_primers[name] = found_pair
            other_options[name] = other_pairs

        else:  # coverage
            pair_candidates_sorted = sorted(
                pair_candidates,
                key=lambda p: (-len(p["coverage"]), p["tm_diff"], p["cross_dimer"], p["prod_dist"])
            )
            best = pair_candidates_sorted[0]
            best_pair = {
                "forward": best["forward"],
                "reverse": best["reverse"],
                "product_size": best["product_size"]
            }
            other_pairs = []
            for alt in pair_candidates_sorted[:10]:
                other_pairs.append({
                    "forward": alt["forward"],
                    "reverse": alt["reverse"],
                    "product_size": alt["product_size"]
                })
            print(f"     ✅ Best pair covers {len(best['coverage'])} isoforms, Tm diff={best['tm_diff']:.1f}°C")
            best_primers[name] = best_pair
            other_options[name] = other_pairs

    return {"best": best_primers, "other": other_options}

# ============================================================================== 
# RACE design (unchanged) 
# ============================================================================== 
def design_race_primers(sequences, params, mode):
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    MAX_RELAX = 4
    tm_mid = (params["tm_min"] + params["tm_max"]) / 2.0
    valid_sequences = [(rec.id, str(rec.seq)) for rec in sequences
                       if len(str(rec.seq)) >= 2 * params["window_size"]]
    if not valid_sequences:
        return {"best": best_primers, "other": other_options}
    for name, seq in valid_sequences:
        seq = seq.upper()
        fw_region = seq[:params["window_size"]]
        rev_region = seq[-params["window_size"]:]
        fw_candidates = find_candidate_primers(
            fw_region, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )
        rev_template_candidates = find_candidate_primers(
            rev_region, params["k_range"], params["gc_min"], params["gc_max"],
            params["tm_min"], params["tm_max"]
        )
        if not fw_candidates or not rev_template_candidates:
            best_primers[name] = None
            other_options[name] = []
            continue
        fw_candidates = list(dict.fromkeys([fw.upper() for fw in fw_candidates]))
        rev_candidates = list(dict.fromkeys([rv.upper() for rv in rev_template_candidates]))
        fw_candidates_sorted = sorted(
            fw_candidates,
            key=lambda k: abs(calculate_melting_temp(k) - tm_mid)
        )
        rev_candidates_sorted = sorted(
            rev_candidates,
            key=lambda k: abs(calculate_melting_temp(k) - tm_mid)
        )
        print(f"🔍 [{name}] RACE - After Tm pruning: {len(fw_candidates_sorted)} FW, {len(rev_candidates_sorted)} REV")
        if mode == "specificity":
            pair_candidates = []
            for fw in fw_candidates_sorted:
                fw_pos = seq.find(fw)
                if fw_pos == -1:
                    continue
                for rv_template in rev_candidates_sorted:
                    rv_template_pos = seq.find(rv_template)
                    if rv_template_pos == -1:
                        continue
                    if fw_pos >= rv_template_pos:
                        continue
                    product_size = (rv_template_pos + len(rv_template)) - fw_pos
                    coverage = get_pair_coverage_set_race(fw, rv_template, sequences)
                    if name not in coverage:
                        continue
                    off_targets = len(coverage) - 1
                    rv_rc = reverse_complement(rv_template)
                    tm_fw = calculate_melting_temp(fw)
                    tm_rv = calculate_melting_temp(rv_rc)
                    tm_diff = abs(tm_fw - tm_rv)
                    cross_dimer = find_longest_complementary_run(fw, reverse_complement(rv_rc))
                    pair_candidates.append({
                        "forward": fw,
                        "reverse": rv_rc,
                        "product_size": product_size,
                        "coverage": coverage,
                        "off_targets": off_targets,
                        "tm_diff": tm_diff,
                        "cross_dimer": cross_dimer
                    })
            print(f"     📊 Generated {len(pair_candidates)} valid RACE pair candidates")
            if not pair_candidates:
                best_primers[name] = None
                other_options[name] = []
                continue
            found_pair = None
            other_pairs = []
            for allowed_off in range(0, MAX_RELAX + 1):
                filtered = [p for p in pair_candidates if p["off_targets"] <= allowed_off]
                if not filtered:
                    continue
                print(f"     🎯 Relaxation level {allowed_off}: {len(filtered)} pairs available")
                filtered_sorted = sorted(
                    filtered,
                    key=lambda p: (p["off_targets"], p["tm_diff"], p["cross_dimer"])
                )
                best = filtered_sorted[0]
                found_pair = {
                    "forward": best["forward"],
                    "reverse": best["reverse"],
                    "product_size": best["product_size"]
                }
                print(f"     ✅ Best pair: {best['off_targets']} off-targets, Tm diff={best['tm_diff']:.1f}°C")
                for alt in filtered_sorted[:10]:
                    other_pairs.append({
                        "forward": alt["forward"],
                        "reverse": alt["reverse"],
                        "product_size": alt["product_size"]
                    })
                break
            best_primers[name] = found_pair
            other_options[name] = other_pairs
        else:
            fw_checks = {fw: check_cached_primer_hits(fw, sequences_tuple) for fw in fw_candidates_sorted}
            rev_checks = {rv: check_cached_primer_hits(rv, sequences_tuple) for rv in rev_candidates_sorted}
            valid_fw = fw_candidates_sorted
            valid_rev = rev_candidates_sorted
            best_pair = None
            other_pairs = []
            for fw, rv in product(valid_fw, valid_rev):
                primer_pair = {"forward": fw, "reverse": rv, "product_size": None}
                if not best_pair:
                    best_pair = primer_pair
                if len(other_pairs) < 10:
                    other_pairs.append(primer_pair)
            best_primers[name] = best_pair
            other_options[name] = other_pairs
    return {"best": best_primers, "other": other_options}

# ============================================================================== 
# SET COVER (unchanged) 
# ============================================================================== 
def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    all_isoform_names = {rec.id for rec in sequences}
    params = params_qpcr if primer_type == "qPCR" else params_race
    def get_regions(seq_str):
        if primer_type == "qPCR":
            return seq_str, seq_str
        else:
            return seq_str[:params["window_size"]], seq_str[-params["window_size"]:]
    fw_kmer_hits = {}
    rv_kmer_hits = {}
    for rec in sequences:
        seq_str = str(rec.seq)
        fw_region, rv_region = get_regions(seq_str)
        for kmer in extract_kmers(fw_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"]:
                fw_kmer_hits.setdefault(kmer, set()).add(rec.id)
        for kmer in extract_kmers(rv_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"]:
                rv_kmer_hits.setdefault(kmer, set()).add(rec.id)
    selected_pairs = []
    pending_isoforms = all_isoform_names.copy()
    combinations = []
    for fw_kmer, fw_hits in fw_kmer_hits.items():
        for rv_kmer, rv_hits in rv_kmer_hits.items():
            common_hits = fw_hits.intersection(rv_hits)
            if common_hits:
                combinations.append((fw_kmer, rv_kmer, common_hits))
    while pending_isoforms and combinations:
        best_pair = max(combinations, key=lambda item: len(item[2].intersection(pending_isoforms)))
        fw, rv, covered_by_best = best_pair
        newly_covered = covered_by_best.intersection(pending_isoforms)
        if not newly_covered:
            break
        selected_pairs.append({
            "forward": fw,
            "reverse": rv,
            "isoforms": list(newly_covered)
        })
        pending_isoforms -= newly_covered
        combinations.remove(best_pair)
    return {"selected": selected_pairs, "pending": list(pending_isoforms)}

# ============================================================================== 
# SCRIPT / IO HELPERS 
# ============================================================================== 
def parse_params(qpcr_params, race_params):
    if qpcr_params:
        qpcr_params = {
            "k_range": range(qpcr_params["k_range"][0], qpcr_params["k_range"][1]),   
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
            "k_range": range(race_params["k_range"][0], race_params["k_range"][1]),   
            "gc_min": race_params["gc_range"][0]*100,
            "gc_max": race_params["gc_range"][1]*100,
            "tm_min": race_params["tm_range"][0],
            "tm_max": race_params["tm_range"][1],
        }
    return qpcr_params, race_params

def create_design(fasta_file_path, primer_type, design_mode, qpcr_params, race_params):
    print(f"📂 Reading FASTA file: {fasta_file_path}")
    num_sequences = 0
    try:
        sequence_set = list(SeqIO.parse(fasta_file_path, "fasta"))
        if not sequence_set:
            print("❌ Error: No sequences found in the FASTA file.")
            return (None, None, None)
        num_sequences = len(sequence_set)
        print(f"✅ Loaded {num_sequences} sequences.")
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at '{fasta_file_path}'")
        return (None, None, None)
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
        for item in results["selected"]:
            coverage_df_rows.append({
                "Primer_Forward": item["forward"],
                "Primer_Reverse": item["reverse"],
                "Covered_Isoforms": ";".join(item["isoforms"]),
            })
        coverage_df = pd.DataFrame(coverage_df_rows) if coverage_df_rows else None
        if coverage_df is not None and not coverage_df.empty:
            coverage_df.index = range(1, len(coverage_df) + 1)
        print("✅ Finished coverage design.")
        return num_sequences, coverage_df, None

if __name__ == "__main__":
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
