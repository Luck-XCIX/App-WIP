# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import logging
import re
from functools import lru_cache
from itertools import product
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# ---------------------------------------------------------------------------
# Third-party
# ---------------------------------------------------------------------------
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
__all__ = [
    "create_design",
    "design_qpcr_primers",
    "design_race_primers",
    "find_coverage_primers",
    "parse_params",
]

# ---------------------------------------------------------------------------
# Configure module logger
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(handler)
logger.setLevel(logging.INFO)

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
COMPLEMENT_MAP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def calculate_gc_content(sequence: str) -> float:
    """Return GC percent of `sequence` (0.0-100.0), rounded to 2 decimals."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc = seq.count("G") + seq.count("C")
    return round(gc / len(seq) * 100.0, 2)


def calculate_melting_temp(sequence: str, na_conc: float = 50e-3) -> float:
    """
    Calculate Tm using Biopython nearest-neighbor (Tm_NN).
    Returns 0.0 on failure.
    """
    try:
        return round(mt.Tm_NN(str(Seq(sequence)), Na=na_conc), 1)
    except Exception:
        return 0.0


def extract_kmers(sequence: str, k_range: range) -> Iterable[str]:
    """Yield all k-mers for each k in `k_range` from `sequence`."""
    if not sequence:
        return
    seq_len = len(sequence)
    for k in k_range:
        if 0 < k <= seq_len:
            for i in range(seq_len - k + 1):
                yield sequence[i : i + k]


def reverse_complement(sequence: str) -> str:
    """Return reverse complement using Biopython Seq utilities."""
    return str(Seq(sequence).reverse_complement())


def are_bases_complementary(base1: str, base2: str) -> bool:
    """Fast complement check using lookup table."""
    return COMPLEMENT_MAP.get(base1.upper()) == base2.upper()


def find_longest_complementary_run(seq1: str, seq2: str) -> int:
    """
    Return the length of the longest perfectly complementary substring between two sequences
    (similar to longest common substring but using complementarity).
    Uses a dynamic programming approach optimized for memory.
    """
    m, n = len(seq1), len(seq2)
    if m == 0 or n == 0:
        return 0

    dp_row = [0] * (n + 1)
    max_len = 0
    for i in range(1, m + 1):
        prev_diagonal = 0
        for j in range(1, n + 1):
            cur = dp_row[j]
            if are_bases_complementary(seq1[i - 1], seq2[j - 1]):
                dp_row[j] = prev_diagonal + 1
                if dp_row[j] > max_len:
                    max_len = dp_row[j]
            else:
                dp_row[j] = 0
            prev_diagonal = cur
    return max_len


def extract_isoform_id(name: str) -> Optional[str]:
    """
    Extract isoform ID from a sequence name.
    Example: 'CASC15-205 cDNA (4297 bp)' -> '205'
    Returns None if no match is found.
    """
    match = re.search(r"[A-Za-z0-9_-]+-(\d+)", name)
    if not match:
        logger.debug("Could not extract isoform ID from: %s", name)
        return None
    return match.group(1)


# ---------------------------------------------------------------------------
# Primer hit checking (cached)
# ---------------------------------------------------------------------------
def check_primer_hits(primer: str, sequences: Sequence[SeqRecord]) -> Dict[str, object]:
    """
    Return a dict describing where `primer` (or its reverse complement) is found.
    Returns:
      {
        "total": int,
        "isoforms": List[str]
      }
    """
    primer_rc = reverse_complement(primer)
    hits: List[str] = []
    for rec in sequences:
        seq_str = str(rec.seq)
        if primer in seq_str or primer_rc in seq_str:
            hits.append(rec.id)
    return {"total": len(hits), "isoforms": hits}


@lru_cache(maxsize=10_000)
def check_cached_primer_hits(primer: str, sequences_tuple: Tuple[Tuple[str, str], ...]) -> Dict[str, object]:
    """
    Cached wrapper for check_primer_hits. Accepts `sequences_tuple` which should be:
    tuple((id, seq_str), ...)
    """
    # Convert back to SeqRecord objects
    records = [SeqRecord(Seq(seq), id=rid) for rid, seq in sequences_tuple]
    return check_primer_hits(primer, records)


def prepare_sequences_for_caching(sequences: Sequence[SeqRecord]) -> Tuple[Tuple[str, str], ...]:
    """Convert SeqRecord list to a hashable tuple suitable for caching."""
    return tuple((rec.id, str(rec.seq)) for rec in sequences)


# ---------------------------------------------------------------------------
# Candidate primer discovery
# ---------------------------------------------------------------------------
def find_candidate_primers(sequence: str, k_range: range, gc_min: float, gc_max: float, tm_min: float, tm_max: float) -> List[str]:
    """
    Return candidate k-mers from `sequence` that satisfy GC and Tm constraints.
    GC is expected as percentage (0-100).
    """
    candidates: List[str] = []
    seq = sequence.upper()
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

# ---------------------------------------------------------------------------
# Design functions
# ---------------------------------------------------------------------------
def design_race_primers(sequences: Sequence[SeqRecord], params: Dict, mode: str) -> Dict[str, List[Dict]]:
    """
    Design primers for RACE mode.
    - sequences: list of SeqRecord
    - params: dictionary with keys window_size, k_range (range), gc_min, gc_max, tm_min, tm_max
    - mode: "specificity" or "coverage" (this function currently supports specificity-like behavior)
    Returns: {"best": {isoform: pair or None}, "other": {isoform: [pairs]}}
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers: Dict[str, Optional[Dict]] = {}
    other_options: Dict[str, List[Dict]] = {}

    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()

        if len(seq) < 2 * params["window_size"]:
            best_primers[name] = None
            other_options[name] = []
            continue

        fw_region = seq[: params["window_size"]]
        rev_region = seq[-params["window_size"] :]
        rev_template = reverse_complement(rev_region)

        fw_candidates = find_candidate_primers(
            fw_region, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"]
        )
        rev_candidates = find_candidate_primers(
            rev_template, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"]
        )

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
            if isoform_id:
                valid_fw = [fw for fw, chk in fw_checks.items() if isoform_id in chk["isoforms"]]
                valid_rev = [rv for rv, chk in rev_checks.items() if isoform_id in chk["isoforms"]]

        best_pair = None
        other_pairs: List[Dict] = []

        for fw, rv in product(valid_fw, valid_rev):
            pair = {"forward": fw, "reverse": rv, "product_size": None}
            if mode == "specificity":
                if not best_pair:
                    best_pair = pair
                if len(other_pairs) < 10:
                    other_pairs.append(pair)

        best_primers[name] = best_pair
        other_options[name] = other_pairs

    return {"best": best_primers, "other": other_options}


def design_qpcr_primers(sequences: Sequence[SeqRecord], params: Dict, mode: str) -> Dict[str, List[Dict]]:
    """
    Design qPCR primers by searching across the entire sequence.
    - params expects: k_range (range), gc_min, gc_max, tm_min, tm_max, prod_min, prod_max
    Returns: {"best": {...}, "other": {...}}
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers: Dict[str, Optional[Dict]] = {}
    other_options: Dict[str, List[Dict]] = {}

    for rec in sequences:
        name = rec.id
        seq = str(rec.seq).upper()

        if len(seq) < params["prod_min"]:
            best_primers[name] = None
            other_options[name] = []
            continue

        # generate forward and reverse template candidates from full sequence
        fw_candidates = find_candidate_primers(seq, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])
        rev_template_candidates = find_candidate_primers(seq, params["k_range"], params["gc_min"], params["gc_max"], params["tm_min"], params["tm_max"])

        # mapping from actual reverse primer (RC) -> template k-mer found in sequence
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
            if isoform_id:
                valid_fw = [fw for fw, chk in fw_checks.items() if isoform_id in chk["isoforms"]]
                valid_rev = [rv for rv, chk in rev_checks.items() if isoform_id in chk["isoforms"]]

        logger.debug("[%s] FW candidates=%d REV candidates=%d", name, len(valid_fw), len(valid_rev))

        best_pair = None
        other_pairs: List[Dict] = []

        for fw, rv in product(valid_fw, valid_rev):
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

            pair = {"forward": fw, "reverse": rv, "product_size": product_size}
            if mode == "specificity":
                if not best_pair:
                    best_pair = pair
                if len(other_pairs) < 10:
                    other_pairs.append(pair)

        best_primers[name] = best_pair
        other_options[name] = other_pairs

    return {"best": best_primers, "other": other_options}


def find_coverage_primers(sequences: Sequence[SeqRecord], params_qpcr: Dict, params_race: Dict, primer_type: str) -> Dict[str, object]:
    """
    Greedy set-cover selection of primer pairs that together amplify (cover) as many isoforms as possible.

    Returns:
      {
        "selected": [ {"forward": f, "reverse": r, "isoforms": [id,...]}, ... ],
        "pending": [isoform_id, ...]
      }
    """
    all_isoform_names: Set[str] = {rec.id for rec in sequences}
    params = params_qpcr if primer_type.lower() == "qpcr" else params_race

    # Compute valid k-mers per isoform
    fw_kmers_by_isoform: Dict[str, Set[str]] = {}
    rv_kmers_by_isoform: Dict[str, Set[str]] = {}

    for rec in sequences:
        seq_str = str(rec.seq).upper()
        if primer_type.lower() == "qpcr":
            fw_region = seq_str
            rv_region = seq_str
        else:
            fw_region = seq_str[: params["window_size"]]
            rv_region = seq_str[-params["window_size"] :]

        fw_valid: Set[str] = set()
        rv_valid: Set[str] = set()

        for kmer in extract_kmers(fw_region, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(kmer) <= params["tm_max"]:
                fw_valid.add(kmer)

        rv_region_to_use = rv_region if primer_type.lower() == "qpcr" else reverse_complement(rv_region)
        for kmer in extract_kmers(rv_region_to_use, params["k_range"]):
            if params["gc_min"] <= calculate_gc_content(kmer) <= params["gc_max"] and params["tm_min"] <= calculate_melting_temp(kmer) <= params["tm_max"]:
                rv_valid.add(kmer)

        fw_kmers_by_isoform[rec.id] = fw_valid
        rv_kmers_by_isoform[rec.id] = rv_valid

    all_fw_primers = set().union(*fw_kmers_by_isoform.values()) if fw_kmers_by_isoform else set()
    all_rv_primers = set().union(*rv_kmers_by_isoform.values()) if rv_kmers_by_isoform else set()

    logger.info("Total unique FW primers: %d", len(all_fw_primers))
    logger.info("Total unique RV primers: %d", len(all_rv_primers))

    def get_pair_coverage(fw_primer: str, rv_primer: str) -> Set[str]:
        covered: Set[str] = set()
        for rec in sequences:
            isoform_id = rec.id
            seq_str = str(rec.seq).upper()

            if fw_primer not in fw_kmers_by_isoform.get(isoform_id, set()) or rv_primer not in rv_kmers_by_isoform.get(isoform_id, set()):
                continue

            fw_pos = seq_str.find(fw_primer)

            if primer_type.lower() == "qpcr":
                rv_rc = reverse_complement(rv_primer)
                rv_pos = seq_str.find(rv_rc)
                if fw_pos != -1 and rv_pos != -1 and fw_pos < rv_pos:
                    product_size = (rv_pos + len(rv_rc)) - fw_pos
                    if params["prod_min"] <= product_size <= params["prod_max"]:
                        covered.add(isoform_id)
            else:
                # RACE: ensure forward exists and reverse exists in 3' window
                if fw_pos != -1:
                    region_3prime = seq_str[-params["window_size"] :]
                    rv_rc = reverse_complement(rv_primer)
                    if rv_rc in region_3prime:
                        covered.add(isoform_id)
        return covered

    # Build candidate list with coverage sets
    combinations: List[Tuple[str, str, Set[str]]] = []
    for fw in all_fw_primers:
        for rv in all_rv_primers:
            cov = get_pair_coverage(fw, rv)
            if cov:
                combinations.append((fw, rv, cov))

    logger.info("Valid primer pair combinations: %d", len(combinations))

    # Greedy set cover
    selected_pairs: List[Dict] = []
    pending_isoforms = set(all_isoform_names)

    while pending_isoforms and combinations:
        best_pair = max(combinations, key=lambda item: len(item[2].intersection(pending_isoforms)))
        fw, rv, covered_by_best = best_pair
        newly_covered = covered_by_best.intersection(pending_isoforms)
        if not newly_covered:
            break
        selected_pairs.append({"forward": fw, "reverse": rv, "isoforms": sorted(list(covered_by_best))})
        logger.info("Selected pair covers %d new isoforms (total %d)", len(newly_covered), len(covered_by_best))
        pending_isoforms -= newly_covered
        combinations.remove(best_pair)

    if pending_isoforms:
        logger.warning("Could not cover %d isoforms: %s", len(pending_isoforms), ", ".join(sorted(pending_isoforms)))

    return {"selected": selected_pairs, "pending": sorted(list(pending_isoforms))}


# ---------------------------------------------------------------------------
# Parameter parsing and top-level pipeline
# ---------------------------------------------------------------------------
def parse_params(qpcr_params: Optional[Dict], race_params: Optional[Dict]) -> Tuple[Optional[Dict], Optional[Dict]]:
    """
    Normalize parameter structures coming from a UI (e.g. ranges given as tuples).
    Expected input format (examples):
      qpcr_params = {
          "k_range": (18, 25),
          "gc_range": (0.4, 0.6),
          "tm_range": (55, 65),
          "prod_range": (70, 200)
      }
      race_params = {
          "window_size": 200,
          "k_range": (18, 25),
          "gc_range": (0.4, 0.6),
          "tm_range": (55, 65)
      }
    Returns normalized dicts with percentages and range() objects.
    """
    def _normalize_qpcr(p: Dict) -> Dict:
        return {
            "k_range": range(p["k_range"][0], p["k_range"][1]),
            "gc_min": p["gc_range"][0] * 100,
            "gc_max": p["gc_range"][1] * 100,
            "prod_min": p["prod_range"][0],
            "prod_max": p["prod_range"][1],
            "tm_min": p["tm_range"][0],
            "tm_max": p["tm_range"][1],
        }

    def _normalize_race(p: Dict) -> Dict:
        return {
            "window_size": p["window_size"],
            "k_range": range(p["k_range"][0], p["k_range"][1]),
            "gc_min": p["gc_range"][0] * 100,
            "gc_max": p["gc_range"][1] * 100,
            "tm_min": p["tm_range"][0],
            "tm_max": p["tm_range"][1],
        }

    qpcr_out = _normalize_qpcr(qpcr_params) if qpcr_params else None
    race_out = _normalize_race(race_params) if race_params else None
    return qpcr_out, race_out


def create_design(
    fasta_file_path: str,
    primer_type: str,
    design_mode: str,
    qpcr_params: Optional[Dict],
    race_params: Optional[Dict],
) -> Tuple[Optional[int], Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Top-level pipeline:
      - fasta_file_path: path to FASTA
      - primer_type: "qPCR" or "RACE"
      - design_mode: "specificity" or "coverage"
      - qpcr_params, race_params: parameter dicts as described in parse_params docstring

    Returns:
      (num_sequences, results_df_1, results_df_2)

    For specificity:
      returns (num_sequences, best_df, other_df)
    For coverage:
      returns (num_sequences, coverage_df, None)
    """
    logger.info("Reading FASTA: %s", fasta_file_path)
    try:
        records = list(SeqIO.parse(fasta_file_path, "fasta"))
    except FileNotFoundError:
        logger.error("FASTA file not found: %s", fasta_file_path)
        return None, None, None

    if not records:
        logger.error("No sequences found in FASTA: %s", fasta_file_path)
        return None, None, None

    num_sequences = len(records)
    logger.info("Loaded %d sequences.", num_sequences)

    qpcr_params_norm, race_params_norm = parse_params(qpcr_params, race_params)

    if design_mode.lower() == "specificity":
        logger.info("Running SPECIFICITY design for %s primers", primer_type)
        if primer_type.lower() == "qpcr":
            results = design_qpcr_primers(records, qpcr_params_norm, design_mode)
        elif primer_type.lower() == "race":
            results = design_race_primers(records, race_params_norm, design_mode)
        else:
            logger.error("Unsupported primer_type: %s", primer_type)
            return None, None, None

        # Build result DataFrames
        best_rows: List[Dict] = []
        other_rows: List[Dict] = []

        for isoform, pair in results["best"].items():
            if not pair:
                continue
            f, r = pair["forward"], pair["reverse"]
            best_rows.append(
                {
                    "Isoform": isoform,
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
                }
            )

        for isoform, pairs in results["other"].items():
            for pair in pairs:
                if not pair:
                    continue
                f, r = pair["forward"], pair["reverse"]
                other_rows.append(
                    {
                        "Isoform": isoform,
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
                    }
                )

        best_df = pd.DataFrame(best_rows) if best_rows else None
        other_df = pd.DataFrame(other_rows) if other_rows else None

        if best_df is not None and not best_df.empty:
            best_df.index = range(1, len(best_df) + 1)
        if other_df is not None and not other_df.empty:
            other_df.index = range(1, len(other_df) + 1)

        logger.info("Specificity design complete.")
        return num_sequences, best_df, other_df

    elif design_mode.lower() == "coverage":
        logger.info("Running COVERAGE design for %s primers", primer_type)
        results = find_coverage_primers(records, qpcr_params_norm, race_params_norm, primer_type)
        coverage_rows: List[Dict] = []
        for idx, item in enumerate(results["selected"], 1):
            f, r = item["forward"], item["reverse"]
            coverage_rows.append(
                {
                    "Pair_ID": f"Pair_{idx}",
                    "Primer_Forward": f,
                    "Primer_Reverse": r,
                    "GC_Forward": calculate_gc_content(f),
                    "GC_Reverse": calculate_gc_content(r),
                    "Tm_Forward": calculate_melting_temp(f),
                    "Tm_Reverse": calculate_melting_temp(r),
                    "Num_Isoforms": len(item["isoforms"]),
                    "Covered_Isoforms": ";".join(item["isoforms"]),
                }
            )

        coverage_df = pd.DataFrame(coverage_rows) if coverage_rows else None
        if coverage_df is not None and not coverage_df.empty:
            coverage_df.index = range(1, len(coverage_df) + 1)

        logger.info("Coverage design complete: %d pairs found, %d pending", len(results["selected"]), len(results["pending"]))
        return num_sequences, coverage_df, None

    else:
        logger.error("Unsupported design_mode: %s", design_mode)
        return None, None, None


# ---------------------------------------------------------------------------
# Example usage when run as a script (very small demo)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    demo_fasta = "example.fasta"  # replace with a real FASTA for local testing

    # Example parameter shapes expected by parse_params
    qpcr_example = {
        "k_range": (18, 21),
        "gc_range": (0.4, 0.6),
        "tm_range": (55, 65),
        "prod_range": (70, 200),
    }
    race_example = {
        "window_size": 200,
        "k_range": (18, 21),
        "gc_range": (0.4, 0.6),
        "tm_range": (55, 65),
    }

    # Demo: change file path and parameters to run locally
    ns, best_df, other_df = create_design(demo_fasta, primer_type="qPCR", design_mode="specificity", qpcr_params=qpcr_example, race_params=race_example)
    logger.info("Result: %s sequences processed", ns)
    if best_df is not None:
        logger.info("Best primers:\n%s", best_df.head().to_string())
