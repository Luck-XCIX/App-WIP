# ==============================================================================
# --- CORRECTED PRIMER DESIGN FUNCTIONS ---
# ==============================================================================

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
