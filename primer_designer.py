# ==============================================================================
# --- CORRECTED PRIMER DESIGN FUNCTIONS ---
# ==============================================================================

def design_race_primers(sequences, params, mode):
    """
    Main logic for designing RACE primers for either specificity or coverage.
    """
    sequences_tuple = prepare_sequences_for_caching(sequences)
    best_primers = {}
    other_options = {}
    
    for rec in sequences:
        name = rec.id
        seq = str(rec.seq)
        
        if len(seq) < 2 * params["window_size"]:
            best_primers[name] = None
            other_options[name] = []
            continue
        
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


def design_qpcr_primers(sequences, params, mode):
    """
    Main logic for designing qPCR primers for either specificity or coverage.
    For qPCR: searches primers across THE ENTIRE SEQUENCE, not limited regions.
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
        
        # CORRECCIÓN: Buscar en TODA la secuencia, no en regiones limitadas
        fw_candidates = find_candidate_primers(
            seq,  # Toda la secuencia
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        
        rev_template_candidates = find_candidate_primers(
            seq,  # Toda la secuencia
            params["k_range"],
            params["gc_min"],
            params["gc_max"],
            params["tm_min"],
            params["tm_max"]
        )
        
        # Map reverse complement for actual primers
        rev_candidates_map = {reverse_complement(p): p for p in rev_template_candidates}
        
        if not fw_candidates or not rev_candidates_map:
            best_primers[name] = None
            other_options[name] = []
            continue
        
        # Check which isoforms each primer hits
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


def parse_params(qpcr_params, race_params):
    """Parse and format parameters for primer design."""
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
            return (None, None, None) if design_mode == "specificity" else (None, None, None)
        num_sequences = len(sequence_set)
        print(f"✅ Loaded {num_sequences} sequences.")
    except FileNotFoundError:
        print(f"❌ Error: FASTA file not found at '{fasta_file_path}'")
        return (None, None, None) if design_mode == "specificity" else (None, None, None)
    
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
        return (None, None, None) if design_mode == "specificity" else (None, None, None)
    
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


def find_coverage_primers(sequences, params_qpcr, params_race, primer_type):
    """
    Finds a minimal set of primer pairs to amplify all isoforms (Set Cover Problem).
    """
    all_isoform_names = {rec.id for rec in sequences}
    params = params_qpcr if primer_type == "qPCR" else params_race
    
    # Pre-compute valid k-mers per isoform
    fw_kmers_by_isoform = {}
    rv_kmers_by_isoform = {}
    
    for rec in sequences:
        seq_str = str(rec.seq)
        
        if primer_type == "qPCR":
            # CORRECCIÓN: Para qPCR, buscar en toda la secuencia
            fw_region = seq_str
            rv_region = seq_str
        else:  # RACE
            # Para RACE, usar ventanas específicas
            fw_region = seq_str[:params["window_size"]]
            rv_region = seq_str[-params["window_size"]:]
        
        fw_valid = set()
        rv_valid = set()
        
        # Extract valid forward k-mers
        for kmer in extract_kmers(fw_region, params["k_range"]):
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if (params["gc_min"] <= gc <= params["gc_max"] and 
                params["tm_min"] <= tm <= params["tm_max"]):
                fw_valid.add(kmer)
        
        # Extract valid reverse k-mers
        if primer_type == "qPCR":
            # Para qPCR: buscar en toda la secuencia (template strand)
            rv_region_to_use = rv_region
        else:  # RACE
            # Para RACE: reverse complement de la ventana 3'
            rv_region_to_use = reverse_complement(rv_region)
        
        for kmer in extract_kmers(rv_region_to_use, params["k_range"]):
            gc = calculate_gc_content(kmer)
            tm = calculate_melting_temp(kmer)
            if (params["gc_min"] <= gc <= params["gc_max"] and 
                params["tm_min"] <= tm <= params["tm_max"]):
                rv_valid.add(kmer)
        
        fw_kmers_by_isoform[rec.id] = fw_valid
        rv_kmers_by_isoform[rec.id] = rv_valid
    
    # Build all valid primer pairs and their coverage
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
                    # Check reverse exists in 3' region
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
