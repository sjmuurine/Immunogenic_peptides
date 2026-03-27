import pandas as pd
import csv

def map_epitopes_to_proteins(output_dir):
    """
    Maps epitopes back to their source proteins with full details.
    """
    
    print("="*60)
    print("MAPPING EPITOPES TO SOURCE PROTEINS")
    print("="*60)
    
    # Load the prioritized epitopes (or you can use TOP50, TIER1, etc.)
    print("\nLoading epitope lists...")
    
    # We'll process multiple files
    files_to_process = {
        'TOP50': f"{output_dir}/TOP50_epitopes.txt",
        'TIER1': f"{output_dir}/TIER1_top_priority_epitopes.txt",
        'TIER2': f"{output_dir}/TIER2_high_priority_epitopes.txt"
    }
    
    # Load source data
    print("Loading source data...")
    
    # IEDB data
    iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    
    # K-mer exact data
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    exact_df = pd.read_csv(exact_file, sep='\t')
    
    # K-mer fuzzy detailed data
    fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
    fuzzy_df = pd.read_csv(fuzzy_file)
    
    # Process each epitope list
    for list_name, file_path in files_to_process.items():
        print(f"\n{'='*60}")
        print(f"Processing {list_name}...")
        print(f"{'='*60}")
        
        epitope_list = pd.read_csv(file_path, sep='\t')
        print(f"Epitopes to map: {len(epitope_list)}")
        
        # Create detailed mapping
        all_mappings = []
        
        for idx, row in epitope_list.iterrows():
            seq = row['Sequence']
            epitope_id = row['Epitope_ID']
            source_methods = row['Source_Methods']
            
            proteins_found = []
            
            # 1. Check IEDB matches
            if 'IEDB' in source_methods or 'HIGH_CONFIDENCE' in source_methods:
                iedb_matches = iedb_df[iedb_df['Epitope_Sequence'] == seq]
                
                for _, match in iedb_matches.iterrows():
                    protein_header = match['Protein_Full_Header']
                    data_source = match['Data_Source']
                    
                    # Parse protein header to extract components
                    # Format: ProteinID Species_OG_Conservation_Localization
                    parts = protein_header.split(' ', 1)
                    protein_id = parts[0] if len(parts) > 0 else protein_header
                    
                    if len(parts) > 1:
                        details = parts[1].split('_')
                        species = details[0] if len(details) > 0 else 'Unknown'
                        orthogroup = details[1] if len(details) > 1 else 'Unknown'
                        conservation = details[2] if len(details) > 2 else 'Unknown'
                        localization = details[3] if len(details) > 3 else 'Unknown'
                    else:
                        species = orthogroup = conservation = localization = 'Unknown'
                    
                    proteins_found.append({
                        'Epitope_ID': epitope_id,
                        'Sequence': seq,
                        'Method': 'IEDB',
                        'Data_Source': data_source,
                        'Protein_ID': protein_id,
                        'Species': species,
                        'Orthogroup': orthogroup,
                        'Conservation': conservation,
                        'Localization': localization,
                        'Percent_Identity': match['Percent_Identity'],
                        'Full_Protein_Header': protein_header
                    })
            
            # 2. Check K-mer exact matches
            if 'Kmer_EXACT' in source_methods or 'HIGH_CONFIDENCE' in source_methods:
                exact_matches = exact_df[exact_df['Motif'] == seq]
                
                for _, match in exact_matches.iterrows():
                    seq_ids_str = str(match['Sequence_IDs'])
                    protein_list = [p.strip() for p in seq_ids_str.split(',')]
                    
                    for protein_header in protein_list:
                        # Determine if ANATOMIA or LITERATURE
                        if protein_header.startswith('LITERATURE'):
                            data_source = 'LITERATURE'
                            # Format: LITERATURE_filename_proteinID
                            parts = protein_header.split('_', 2)
                            protein_id = parts[2] if len(parts) > 2 else protein_header
                            species = 'LITERATURE'
                            orthogroup = 'N/A'
                            conservation = 'N/A'
                            localization = 'N/A'
                        else:
                            data_source = 'ANATOMIA'
                            # Format: ProteinID Species_OG_Conservation_Localization
                            parts = protein_header.split(' ', 1)
                            protein_id = parts[0] if len(parts) > 0 else protein_header
                            
                            if len(parts) > 1:
                                details = parts[1].split('_')
                                species = details[0] if len(details) > 0 else 'Unknown'
                                orthogroup = details[1] if len(details) > 1 else 'Unknown'
                                conservation = details[2] if len(details) > 2 else 'Unknown'
                                localization = details[3] if len(details) > 3 else 'Unknown'
                            else:
                                species = orthogroup = conservation = localization = 'Unknown'
                        
                        proteins_found.append({
                            'Epitope_ID': epitope_id,
                            'Sequence': seq,
                            'Method': 'Kmer_EXACT',
                            'Data_Source': data_source,
                            'Protein_ID': protein_id,
                            'Species': species,
                            'Orthogroup': orthogroup,
                            'Conservation': conservation,
                            'Localization': localization,
                            'Percent_Identity': 100,  # Exact match
                            'Full_Protein_Header': protein_header
                        })
            
            # 3. Check K-mer fuzzy matches
            if 'Kmer_FUZZY' in source_methods:
                fuzzy_matches = fuzzy_df[fuzzy_df['Consensus'] == seq]
                
                for _, match in fuzzy_matches.iterrows():
                    protein_header = str(match['Source_Protein'])
                    
                    # Determine if ANATOMIA or LITERATURE
                    if protein_header.startswith('LITERATURE'):
                        data_source = 'LITERATURE'
                        parts = protein_header.split('_', 2)
                        protein_id = parts[2] if len(parts) > 2 else protein_header
                        species = 'LITERATURE'
                        orthogroup = 'N/A'
                        conservation = 'N/A'
                        localization = 'N/A'
                    else:
                        data_source = 'ANATOMIA'
                        parts = protein_header.split(' ', 1)
                        protein_id = parts[0] if len(parts) > 0 else protein_header
                        
                        if len(parts) > 1:
                            details = parts[1].split('_')
                            species = details[0] if len(details) > 0 else 'Unknown'
                            orthogroup = details[1] if len(details) > 1 else 'Unknown'
                            conservation = details[2] if len(details) > 2 else 'Unknown'
                            localization = details[3] if len(details) > 3 else 'Unknown'
                        else:
                            species = orthogroup = conservation = localization = 'Unknown'
                    
                    similarity = match['Similarity'] if 'Similarity' in match else 1.0
                    
                    proteins_found.append({
                        'Epitope_ID': epitope_id,
                        'Sequence': seq,
                        'Method': 'Kmer_FUZZY',
                        'Data_Source': data_source,
                        'Protein_ID': protein_id,
                        'Species': species,
                        'Orthogroup': orthogroup,
                        'Conservation': conservation,
                        'Localization': localization,
                        'Percent_Identity': similarity * 100,
                        'Full_Protein_Header': protein_header
                    })
            
            all_mappings.extend(proteins_found)
        
        # Convert to DataFrame
        mappings_df = pd.DataFrame(all_mappings)
        
        # Save detailed mapping
        output_file = f"{output_dir}/{list_name}_protein_mapping.txt"
        mappings_df.to_csv(output_file, sep='\t', index=False)
        print(f"✓ Saved: {output_file}")
        print(f"  Total protein mappings: {len(mappings_df)}")
        
        # Create summary by epitope
        summary_list = []
        for epitope_id in mappings_df['Epitope_ID'].unique():
            epitope_data = mappings_df[mappings_df['Epitope_ID'] == epitope_id]
            
            anatomia_count = len(epitope_data[epitope_data['Data_Source'] == 'ANATOMIA'])
            literature_count = len(epitope_data[epitope_data['Data_Source'] == 'LITERATURE'])
            
            species_list = epitope_data['Species'].unique()
            og_list = epitope_data[epitope_data['Orthogroup'] != 'N/A']['Orthogroup'].unique()
            
            summary_list.append({
                'Epitope_ID': epitope_id,
                'Sequence': epitope_data.iloc[0]['Sequence'],
                'Total_Proteins': len(epitope_data),
                'ANATOMIA_Proteins': anatomia_count,
                'LITERATURE_Proteins': literature_count,
                'Num_Species': len([s for s in species_list if s not in ['Unknown', 'LITERATURE']]),
                'Species_List': ', '.join([s for s in species_list if s not in ['Unknown', 'LITERATURE']])[:200],
                'Num_Orthogroups': len([o for o in og_list if o != 'Unknown']),
                'Orthogroups_List': ', '.join([o for o in og_list if o != 'Unknown'])[:200]
            })
        
        summary_df = pd.DataFrame(summary_list)
        summary_file = f"{output_dir}/{list_name}_protein_summary.txt"
        summary_df.to_csv(summary_file, sep='\t', index=False)
        print(f"✓ Saved: {summary_file}")
    
    print("\n" + "="*60)
    print("MAPPING COMPLETE")
    print("="*60)
    print("\nGenerated files for each tier:")
    print("  - *_protein_mapping.txt : Detailed protein-level data")
    print("  - *_protein_summary.txt : Summary by epitope")
    print("="*60)

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
map_epitopes_to_proteins(output_dir)