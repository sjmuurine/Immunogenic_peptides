import pandas as pd
import os

def map_reranked_to_full_metadata(cluster_dir, comparison_dir, output_dir):
    """
    Map re-ranked epitope clusters to full protein metadata.
    """
    
    print("="*70)
    print("MAPPING RERANKED EPITOPES TO FULL METADATA")
    print("="*70)
    
    # Load IEDB data
    print("\nLoading IEDB data for links...")
    iedb_links = {}
    try:
        iedb_blast_file = f"{comparison_dir}/../002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt"
        iedb_df = pd.read_csv(iedb_blast_file, sep='\t')
        for _, row in iedb_df.iterrows():
            seq = row['Epitope_Sequence']
            iedb_id = row['Epitope_ID']
            if seq not in iedb_links:
                iedb_links[seq] = []
            if iedb_id not in iedb_links[seq]:
                iedb_links[seq].append(iedb_id)
        print(f"  Loaded IEDB links for {len(iedb_links)} sequences")
    except Exception as e:
        print(f"  Warning: Could not load IEDB links: {e}")
    
    # Load protein mappings
    protein_mappings = {}
    for tier in ['TOP50', 'TIER1', 'TIER2']:
        mapping_file = f"{comparison_dir}/{tier}_protein_mapping.txt"
        try:
            df = pd.read_csv(mapping_file, sep='\t', usecols=['Epitope_ID', 'Data_Source', 'Full_Protein_Header'])
            
            def parse_header(header):
                """
                Parse from Full_Protein_Header.
                Format: ProteinID Genus_species_Orthogroup_Conservation_Localization
                Key insight: Species has underscore (e.g., V_magna), need to find OG to know where species ends
                """
                if pd.isna(header):
                    return 'Unknown', 'Unknown', 'Unknown', 'Unknown', 'Unknown'
                
                header_str = str(header).strip()
                
                # LITERATURE
                if header_str.startswith('LITERATURE'):
                    parts = header_str.split('_')
                    protein_id = header_str.split()[0] if ' ' in header_str else header_str
                    species = 'LITERATURE'
                    
                    for part in parts:
                        part_lower = part.lower()
                        if 'nucleatum' in part_lower:
                            species = 'F_nucleatum'
                        elif 'fragilis' in part_lower:
                            species = 'B_fragilis'
                        elif 'muciniphila' in part_lower:
                            species = 'A_muciniphila'
                        elif 'prausnitzii' in part_lower:
                            species = 'F_prausnitzii'
                        elif 'reuteri' in part_lower:
                            species = 'L_reuteri'
                        elif 'plantarum' in part_lower:
                            species = 'L_plantarum'
                        elif 'casei' in part_lower:
                            species = 'L_casei'
                        elif 'cloacae' in part_lower:
                            species = 'E_cloacae'
                        elif 'faecalis' in part_lower:
                            species = 'E_faecalis'
                        elif 'salivarius' in part_lower:
                            species = 'S_salivarius'
                        elif 'longum' in part_lower:
                            species = 'B_longum'
                        elif 'gasseri' in part_lower:
                            species = 'L_gasseri'
                        elif 'fermentum' in part_lower:
                            species = 'L_fermentum'
                        elif 'adiacens' in part_lower:
                            species = 'G_adiacens'
                        elif 'freudenreichii' in part_lower:
                            species = 'P_freudenreichii'
                    
                    return protein_id, species, 'N/A', 'N/A', 'LITERATURE'
                
                # ANATOMIA: WP_xxxxx.x Genus_species_OGxxxxxx_Conservation_Localization
                try:
                    parts = header_str.split(' ', 1)
                    if len(parts) < 2:
                        return header_str, 'Unknown', 'Unknown', 'Unknown', 'Unknown'
                    
                    protein_id = parts[0]
                    metadata = parts[1]
                    
                    # Find where OG starts (it's the part that starts with "OG" or is "unique")
                    meta_parts = metadata.split('_')
                    
                    og_index = None
                    for i, part in enumerate(meta_parts):
                        if part.startswith('OG') or part == 'unique':
                            og_index = i
                            break
                    
                    if og_index is None or og_index < 1:
                        return protein_id, 'Unknown', 'Unknown', 'Unknown', 'Unknown'
                    
                    # Everything before OG is species (join with underscore)
                    species = '_'.join(meta_parts[:og_index])
                    
                    # OG is at og_index
                    orthogroup = meta_parts[og_index]
                    
                    # Conservation is right after OG
                    if og_index + 1 < len(meta_parts):
                        conservation = meta_parts[og_index + 1]
                    else:
                        conservation = 'Unknown'
                    
                    # Localization is everything after conservation
                    if og_index + 2 < len(meta_parts):
                        localization = ' '.join(meta_parts[og_index + 2:])
                    else:
                        localization = 'Unknown'
                    
                    return protein_id, species, orthogroup, conservation, localization
                
                except Exception as e:
                    return header_str, 'Unknown', 'Unknown', 'Unknown', 'Unknown'
            
            # Parse
            parsed = df['Full_Protein_Header'].apply(parse_header)
            df['Protein_ID'] = [p[0] for p in parsed]
            df['Species'] = [p[1] for p in parsed]
            df['Orthogroup'] = [p[2] for p in parsed]
            df['Conservation'] = [p[3] for p in parsed]
            df['Localization'] = [p[4] for p in parsed]
            
            protein_mappings[tier] = df
            print(f"Loaded {tier}: {len(df)} entries")
            
        except Exception as e:
            print(f"Error loading {tier}: {e}")
            protein_mappings[tier] = pd.DataFrame()
    
    # Process re-ranked files
    datasets = ['TOP50', 'TIER1', 'TIER2', 'COMBINED_TOP_TIERS']
    
    for dataset in datasets:
        reranked_file = f"{cluster_dir}/{dataset}_RERANKED.txt"
        
        try:
            reranked_df = pd.read_csv(reranked_file, sep='\t')
        except:
            print(f"\nSkipping {dataset}")
            continue
        
        print(f"\n{'='*70}")
        print(f"Processing: {dataset}")
        print(f"{'='*70}")
        
        detailed_file = f"{cluster_dir}/{dataset}_cluster_detailed_v2.txt"
        try:
            detailed_df = pd.read_csv(detailed_file, sep='\t')
        except:
            detailed_df = pd.DataFrame()
        
        epitope_metadata = []
        
        for idx, cluster_row in reranked_df.iterrows():
            cluster_id = cluster_row['Cluster_ID']
            rep_sequence = cluster_row['Representative_Sequence']
            
            iedb_link_list = iedb_links.get(rep_sequence, [])
            iedb_link_str = '; '.join(iedb_link_list[:3]) if iedb_link_list else 'N/A'
            
            if len(detailed_df) > 0:
                cluster_epitopes = detailed_df[detailed_df['Cluster_ID'] == cluster_id]
                epitope_ids = cluster_epitopes['Epitope_ID'].tolist()
            else:
                epitope_ids = [cluster_row['Representative_ID']]
            
            all_proteins = []
            for epitope_id in epitope_ids:
                for tier, mapping_df in protein_mappings.items():
                    if len(mapping_df) > 0:
                        matches = mapping_df[mapping_df['Epitope_ID'] == epitope_id]
                        if len(matches) > 0:
                            all_proteins.append(matches)
            
            if not all_proteins:
                epitope_metadata.append({
                    'Cluster_ID': cluster_id,
                    'Rank': idx + 1,
                    'Representative_Sequence': rep_sequence,
                    'Cluster_Size': cluster_row['Cluster_Size'],
                    'Priority_Score': cluster_row['New_Priority_Score'],
                    'Total_Hits': cluster_row['Total_Hits_AllVariants'],
                    'Methods_Used': cluster_row['Methods_Used'],
                    'IEDB_Links': iedb_link_str,
                    'Species_List': 'Not mapped',
                    'Localization_List': 'Not mapped',
                    'Orthogroup_List': 'Not mapped',
                    'Conservation_List': 'Not mapped',
                    'Num_Surface_Proteins': 0,
                    'Suitable_For_Fishing': 'No'
                })
                continue
            
            combined = pd.concat(all_proteins, ignore_index=True)
            combined = combined.drop_duplicates(subset=['Protein_ID'])
            
            def safe_unique(series, exclude=['Unknown', 'N/A', '']):
                values = series.dropna().astype(str).unique()
                filtered = [v for v in values if v not in exclude]
                return sorted(set(filtered))
            
            species = safe_unique(combined['Species'])
            localizations = safe_unique(combined['Localization'], ['Unknown', 'LITERATURE'])
            orthogroups = safe_unique(combined['Orthogroup'])
            conservation = safe_unique(combined['Conservation'])
            
            def is_surface(loc):
                if pd.isna(loc) or str(loc) in ['Unknown', 'LITERATURE']:
                    return False
                loc = str(loc).lower()
                if 'cytoplasmic' in loc and 'membrane' not in loc:
                    return False
                return any(kw in loc for kw in ['membrane', 'extracellular', 'outer', 'surface', 'periplasm', 'secreted', 'cellwall'])
            
            combined['Is_Surface'] = combined['Localization'].apply(is_surface)
            num_surface = combined['Is_Surface'].sum()
            
            anatomia_count = len(combined[combined['Data_Source'] == 'ANATOMIA'])
            literature_count = len(combined[combined['Data_Source'] == 'LITERATURE'])
            
            sample_proteins = combined.head(3)['Full_Protein_Header'].tolist()
            
            epitope_metadata.append({
                'Cluster_ID': cluster_id,
                'Rank': idx + 1,
                'Representative_Sequence': rep_sequence,
                'Sequence_Length': cluster_row['Representative_Length'],
                'Cluster_Size': cluster_row['Cluster_Size'],
                'Priority_Score': cluster_row['New_Priority_Score'],
                'Total_Hits': cluster_row['Total_Hits_AllVariants'],
                'ANATOMIA_Hits': cluster_row['ANATOMIA_Hits_AllVariants'],
                'LITERATURE_Hits': cluster_row['LITERATURE_Hits_AllVariants'],
                'Methods_Used': cluster_row['Methods_Used'],
                'IEDB_Links': iedb_link_str,
                'Num_Species': len(species),
                'Species_List': '; '.join(species) if species else 'None',
                'Num_Localizations': len(localizations),
                'Localization_List': '; '.join(localizations) if localizations else 'Unknown',
                'Num_Orthogroups': len(orthogroups),
                'Orthogroup_List': '; '.join(orthogroups) if orthogroups else 'None',
                'Conservation_List': '; '.join(conservation) if conservation else 'None',
                'Total_Proteins': len(combined),
                'ANATOMIA_Proteins': anatomia_count,
                'LITERATURE_Proteins': literature_count,
                'Num_Surface_Proteins': int(num_surface),
                'Suitable_For_Fishing': 'Yes' if num_surface > 0 else 'No',
                'Human_Identity': cluster_row['Human_Identity'],
                'Sample_Proteins': ' | '.join(sample_proteins)
            })
        
        metadata_df = pd.DataFrame(epitope_metadata)
        
        output_file = f"{output_dir}/{dataset}_FINAL_with_metadata.txt"
        metadata_df.to_csv(output_file, sep='\t', index=False)
        print(f"✓ Saved: {output_file}")
        
        fishable = metadata_df[metadata_df['Suitable_For_Fishing'] == 'Yes']
        if len(fishable) > 0:
            fishable_file = f"{output_dir}/{dataset}_FISHABLE_epitopes.txt"
            fishable.to_csv(fishable_file, sep='\t', index=False)
            print(f"✓ Fishable: {len(fishable)}")
        
        print(f"\nTop 3:")
        for _, row in metadata_df.head(3).iterrows():
            print(f"{row['Rank']}. {row['Representative_Sequence']}")
            print(f"   OG: {row['Orthogroup_List']}")
            print(f"   Conservation: {row['Conservation_List']}")
            print(f"   Localization: {row['Localization_List']}")
    
    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)

# Run
cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
comparison_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'

map_reranked_to_full_metadata(cluster_dir, comparison_dir, output_dir)