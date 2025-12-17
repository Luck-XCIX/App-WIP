# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------

import re
from itertools import product
from functools import lru_cache

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

# -------------------------
# Utility functions
# -------------------------

def calculate_gc_content(seq: str) -> float:
    """Return GC% (0-100) rounded to 2 decimals."""
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return round(gc / len(s) * 100.0, 2)


def extract_kmers(sequence: str, k_range: range):
    """Yield all k-mers for lengths in k_range."""
    if not sequence:
        return
    seq = sequence.upper()
    for k in k_range:
        if 0 < k <= len(seq):
            for i in range(len(seq) - k + 1):
                yield seq[i : i + k]


COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(sequence: str) -> str:
    """Reverse complement using Biopython Seq."""
    return str(Seq(sequence).reverse_complement())


def find_longest_complementary_run(a: str, b: str) -> int:
    """
    Return length of longest perfectly complementary substring between a and b.
    Memory-optimized DP.
    """
    a, b = a.upper(), b.upper()
    m, n = len(a), len(b)
    if m == 0 or n == 0:
        return 0
    dp = [0] * (n + 1)
    best = 0
    for i in range(1, m + 1):
        prev = 0
        for j in range(1, n + 1):
            cur = dp[j]
            if COMPLEMENT_MAP.get(a[i - 1]) == b[j - 1]:
                dp[j] = prev + 1
                if dp[j] > best:
                    best = dp[j]
            else:
                dp[j] = 0
            prev = cur
    return best


def calculate_melting_temp(sequence: str, na_conc: float = 50e-3) -> float:
    """Return Tm using Biopython nearest-neighbor (rounded)."""
    try:
        return round(mt.Tm_NN(str(Seq(sequence)), Na=na_conc), 1)
    except Exception:
        return 0.0


def extract_isoform_id(name: str) -> str:
    """
    Extract isoform numeric ID from a name like 'CASC15-205 ...'.
    Returns empty string if none found.
    """
    m = re.search(r"-([0-9]+)", name)
    return m.group(1) if m else ""


def check_primer_hits(primer: str, sequences):
    """
    Return dict with total hits and comma-separated isoform ids where primer
    or its reverse complement occurs.
    """
    rc = reverse_complement(primer)
    hits = [rec.id for rec in sequences if primer in str(rec.seq).upper() or rc in str(rec.seq).upper()]
    return {"total": len(hits), "isoforms": ",".join(hits)}


def find_candidate_primers(sequence: str, k_range: range, gc_min: float, gc_max: float, tm_min: float, tm_max: float):
    """Return list of k-mers that meet GC% and Tm constraints."""
    seq = sequence.upper()
    candidates = []
    for k in k_range:
        if k <= 0 or k > len(seq):
            continue
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if gc_min <= gc <= gc_max and tm_min <= tm <= tm_max:
                candidates.append(kmer)
    return candidates


@lru_cache(maxsize=10000)
def check_cached_primer_hits(primer: str, sequences_tuple):
    """
    Cached wrapper for check_primer_hits.
    sequences_tuple: tuple of (id, seq_str)
    """
    records = [SeqIO.SeqRecord(Seq(seq), id=name) for name, seq in sequences_tuple]
    return check_primer_hits(primer, records)


def prepare_sequences_for_caching(sequences):
    """Return tuple of (id, seq_str) for caching."""
    return tuple((rec.id, str(rec.seq).upper()) for rec in sequences)


# -------------------------
# Primer design functions
# -------------------------

def design_race_primers(sequences, params, mode):
    """
    Design RACE primers (windowed at ends).
    Returns {"best": {id: pair_or_None}, "other": {id: [pairs]}}.
    """
    seq_tuple = prepare_sequences_for_caching(sequences)
    best = {}
    other = {}

    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()
        if len(seq) < 2 * params["window_size"]:
            best[name] = None
            other[name] = []
            continue

        fw_region = seq[: params["window_size"]]
        rev_region = seq[-params["window_size"] :]
        rev_template = reverse_complement(rev_region)

        fw_cands = find_candidate_primers(fw_region, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])
        rv_cands = find_candidate_primers(rev_template, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])

        if not fw_cands or not rv_cands:
            best[name] = None
            other[name] = []
            continue

        fw_checks = {fw: check_cached_primer_hits(fw, seq_tuple) for fw in fw_cands}
        rv_checks = {rv: check_cached_primer_hits(rv, seq_tuple) for rv in rv_cands}

        valid_fw = fw_cands
        valid_rv = rv_cands
        if mode == "specificity":
            iso = extract_isoform_id(name)
            valid_fw = [fw for fw, chk in fw_checks.items() if iso in chk["isoforms"]]
            valid_rv = [rv for rv, chk in rv_checks.items() if iso in chk["isoforms"]]

        best_pair = None
        other_pairs = []
        for fw, rv in product(valid_fw, valid_rv):
            pair = {"forward": fw, "reverse": rv, "product_size": None}
            if mode == "specificity":
                if not best_pair:
                    best_pair = pair
                if len(other_pairs) < 10:
                    other_pairs.append(pair)

        best[name] = best_pair
        other[name] = other_pairs

    return {"best": best, "other": other}


def design_qpcr_primers(sequences, params, mode):
    """
    Design qPCR primers by scanning full sequences.
    Returns {"best": {id: pair_or_None}, "other": {id: [pairs]}}.
    """
    seq_tuple = prepare_sequences_for_caching(sequences)
    best = {}
    other = {}

    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()
        if len(seq) < params["prod_min"]:
            best[name] = None
            other[name] = []
            continue

        fw_cands = find_candidate_primers(seq, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])
        rev_template_cands = find_candidate_primers(seq, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])

        rev_map = {reverse_complement(p): p for p in rev_template_cands}
        if not fw_cands or not rev_map:
            best[name] = None
            other[name] = []
            continue

        fw_checks = {fw: check_cached_primer_hits(fw, seq_tuple) for fw in fw_cands}
        rev_checks = {rv: check_cached_primer_hits(rv, seq_tuple) for rv in rev_map.keys()}

        valid_fw = fw_cands
        valid_rv = list(rev_map.keys())
        if mode == "specificity":
            iso = extract_isoform_id(name)
            valid_fw = [fw for fw, chk in fw_checks.items() if iso in chk["isoforms"]]
            valid_rv = [rv for rv, chk in rev_checks.items() if iso in chk["isoforms"]]

        best_pair = None
        other_pairs = []
        for fw, rv in product(valid_fw, valid_rv):
            fw_pos = seq.find(fw)
            rv_template = rev_map[rv]
            rv_pos = seq.find(rv_template)
            if fw_pos == -1 or rv_pos == -1 or fw_pos >= rv_pos:
                continue
            product_size = (rv_pos + len(rv_template)) - fw_pos
            if not (params["prod_min"] <= product_size <= params["prod_max"]):
                continue
            pair = {"forward": fw, "reverse": rv, "product_size": product_size}
            if mode == "specificity":
                if not best_pair:
                    best_pair = pair
                if len(other_pairs) < 10:
                    other_pairs.append(pair)

        best[name] = best_pair
        other[name] = other_pairs

    return {"best": best, "other": other}


def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    """
    Greedy set-cover over isoforms: choose primer pairs that cover the most isoforms until none left.
    Returns {"selected": [pairs], "pending": [isoform_ids]}.
    """
    isoforms = {rec.id for rec in sequences}
    params = params_qpcr if primer_type == "qPCR" else params_race

    fw_kmers = {}
    rv_kmers = {}

    for rec in sequences:
        seq = str(rec.seq).upper()
        if primer_type == "qPCR":
            fw_region = rv_region = seq
        else:
            fw_region = seq[: params["window_size"]]
            rv_region = seq[-params["window_size"] :]

        fw_valid = set()
        rv_valid = set()

        for kmer in extract_kmers(fw_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(kmer) <= params["tm_max"]:
                fw_valid.add(kmer)

        rv_region_use = rv_region if primer_type == "qPCR" else reverse_complement(rv_region)
        for kmer in extract_kmers(rv_region_use, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(kmer) <= params["tm_max"]:
                rv_valid.add(kmer)

        fw_kmers[rec.id] = fw_valid
        rv_kmers[rec.id] = rv_valid

    all_fw = set().union(*fw_kmers.values()) if fw_kmers else set()
    all_rv = set().union(*rv_kmers.values()) if rv_kmers else set()

    def pair_coverage(fw, rv):
        covered = set()
        for rec in sequences:
            iso = rec.id
            seq = str(rec.seq).upper()
            if fw not in fw_kmers.get(iso, set()) or rv not in rv_kmers.get(iso, set()):
                continue
            fw_pos = seq.find(fw)
            if primer_type == "qPCR":
                rv_rc = reverse_complement(rv)
                rv_pos = seq.find(rv_rc)
                if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
                    size = (rv_pos + len(rv_rc)) - fw_pos
                    if params["prod_min"] <= size <= params["prod_max"]:
                        covered.add(iso)
            else:
                if fw_pos != -1:
                    region_3 = seq[-params["window_size"] :]
                    if reverse_complement(rv) in region_3:
                        covered.add(iso)
        return covered

    combos = []
    for fw in all_fw:
        for rv in all_rv:
            cov = pair_coverage(fw, rv)
            if cov:
                combos.append((fw, rv, cov))

    selected = []
    pending = set(isoforms)
    while pending and combos:
        best_combo = max(combos, key=lambda x: len(x[2].intersection(pending)))
        fw, rv, cov = best_combo
        new = cov.intersection(pending)
        if not new:
            break
        selected.append({"forward": fw, "reverse": rv, "isoforms": sorted(list(cov))})
        pending -= new
        combos.remove(best_combo)

    return {"selected": selected, "pending": sorted(list(pending))}


# -------------------------
# Helpers & top-level
# -------------------------

def parse_params(qpcr_params, race_params):
    """Normalize UI-style params to internal form (ranges -> range, gc -> percent)."""
    if qpcr_params:
        qpcr_params = {
            "k_range": range(qpcr_params["k_range"][0], qpcr_params["k_range"][1]),
            "gc_min": qpcr_params["gc_range"][0] * 100,
            "gc_max": qpcr_params["gc_range"][1] * 100,
            "prod_min": qpcr_params["prod_range"][0],
            "prod_max": qpcr_params["prod_range"][1],
            "tm_min": qpcr_params["tm_range"][0],
            "tm_max": qpcr_params["tm_range"][1],
        }
    if race_params:
        race_params = {
            "window_size": race_params["window_size"],
            "k_range": range(race_params["k_range"][0], race_params["k_range"][1]),
            "gc_min": race_params["gc_range"][0] * 100,
            "gc_max": race_params["gc_range"][1] * 100,
            "tm_min": race_params["tm_range"][0],
            "tm_max": race_params["tm_range"][1],
        }
    return qpcr_params, race_params


def create_design(fasta_file_path, primer_type, design_mode, qpcr_params, race_params):
    """
    Main pipeline. Returns:
      - specificity: (num_sequences, best_df, other_df)
      - coverage: (num_sequences, coverage_df, None)
    """
    try:
        records = list(SeqIO.parse(fasta_file_path, "fasta"))
    except FileNotFoundError:
        print(f"FASTA not found: {fasta_file_path}")
        return None, None, None

    if not records:
        print("No sequences found in FASTA.")
        return None, None, None

    qpcr_params, race_params = parse_params(qpcr_params, race_params)

    if design_mode == "specificity":
        if primer_type == "qPCR":
            results = design_qpcr_primers(records, qpcr_params, design_mode)
        elif primer_type == "RACE":
            results = design_race_primers(records, race_params, design_mode)
        else:
            print("Unsupported primer type.")
            return None, None, None

        best_rows = []
        other_rows = []
        for iso, pair in results["best"].items():
            if not pair:
                continue
            f, r = pair["forward"], pair["reverse"]
            best_rows.append({
                "Isoform": iso,
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
                "Product_Size": pair.get("product_size"),
            })
        for iso, pairs in results["other"].items():
            for pair in pairs:
                if not pair:
                    continue
                f, r = pair["forward"], pair["reverse"]
                other_rows.append({
                    "Isoform": iso,
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
                    "Product_Size": pair.get("product_size"),
                })

        best_df = pd.DataFrame(best_rows) if best_rows else None
        other_df = pd.DataFrame(other_rows) if other_rows else None
        if best_df is not None and not best_df.empty:
            best_df.index = range(1, len(best_df) + 1)
        if other_df is not None and not other_df.empty:
            other_df.index = range(1, len(other_df) + 1)

        return len(records), best_df, other_df

    elif design_mode == "coverage":
        results = find_coverage_primers(records, qpcr_params, race_params, primer_type)
        rows = []
        for i, item in enumerate(results["selected"], 1):
            f, r = item["forward"], item["reverse"]
            rows.append({
                "Pair_ID": f"Pair_{i}",
                "Primer_Forward": f,
                "Primer_Reverse": r,
                "GC_Forward": calculate_gc_content(f),
                "GC_Reverse": calculate_gc_content(r),
                "Tm_Forward": calculate_melting_temp(f),
                "Tm_Reverse": calculate_melting_temp(r),
                "Num_Isoforms": len(item["isoforms"]),
                "Covered_Isoforms": ";".join(item["isoforms"]),
            })
        cov_df = pd.DataFrame(rows) if rows else None
        if cov_df is not None and not cov_df.empty:
            cov_df.index = range(1, len(cov_df) + 1)
        return len(records), cov_df, None

    else:
        print("Unsupported design mode.")
        return None, None, None
