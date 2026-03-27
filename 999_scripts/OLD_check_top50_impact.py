import pandas as pd
import os

def parse_species_from_protein_id(protein_id):
    """
    Parse species from protein ID - CORRECTED VERSION with LITERATURE support.
    """
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    # Enhanced species patterns
    species_patterns = {
        'fragilis': 'B_fragilis', 'B_fragilis': 'B_fragilis',
        'nucleatum': 'F_nucleatum', 'Fnucleatum': 'F_nucleatum', 'F_nucleatum': 'F_nucleatum',
        'muciniphila': 'A_muciniphila', 'A_muciniphila': 'A_muciniphila',
        'prausnitzii': 'F_prausnitzii', 'F_prausnitzii': 'F_prausnitzii',
        'reuteri': 'L_reuteri', 'L_REUTERI': 'L_reuteri', 'L_reuteri': 'L_reuteri',
        'cloacae': 'E_cloacae', 'Ecloacae': 'E_cloacae', 'E_cloacae': 'E_cloacae',
        'faecalis': 'E_faecalis', 'efaecalis': 'E_faecalis', 'E_faecalis': 'E_faecalis',
        'salivarius': 'S_salivarius', 'S_salivarius': 'S_salivarius',
        'longum': 'B_longum', 'B_longum': 'B_longum', 'B_Longum': 'B_longum',
        'plantarum': 'L_plantarum', 'L_plantarum': 'L_plantarum',
        'casei': 'L_casei', 'Casei': 'L_casei', 'L_casei': 'L_casei',
        'sanguinis': 'T_sanguinis', 'T_sanguinis': 'T_sanguinis',
        'magna': 'V_magna', 'V_magna': 'V_magna',
        'russellii': 'P_russellii', 'P_russellii': 'P_russellii',
        'freudenreichii': 'P_freudenreichii', 'P_freudenreichii': 'P_freudenreichii',
        'fermentum': 'L_fermentum', 'L_fermentum': 'L_fermentum',
        'adiacens': 'G_adiacens', 'G_adiacens': 'G_adiacens',
        'gasseri': 'L_gasseri', 'L_gasseri': 'L_gasseri'
    }
    
    # LITERATURE format
    if 'LITERATURE' in protein_id:
        for pattern, species in species_patterns.items():
            if pattern.lower() in protein_id.lower():
                return species
        return None
    
    # ANATOMIA format
    if ' ' not in protein_id:
        return None
    
    parts = protein_id.split(' ', 1)[1].split('_')
    og_idx = None
    for i, part in enumerate(parts):
        if part.startswith('OG') or part == 'unique':
            og_idx = i
            break
    
    if og_idx and og_idx > 0:
        return '_'.join(parts[:og_idx])
    
    return None

def main():
    """
    Comprehensive comparison: OLD vs NEW species counts for TOP50.
    """
    
    print("="*80)
    print("COMPREHENSIVE TOP50 SPECIES COUNT COMPARISON")
    print("="*80)
    
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    
    # 1. Load cluster-epitope mapping
    detailed_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_cluster_detailed_v2.txt'
    detailed_df = pd.read_csv(detailed_file, sep='\t')
    
    cluster_to_epitopes = {}
    for cluster_id in detailed_df['Cluster_ID'].unique():
        epitope_ids = detailed_df[detailed_df['Cluster_ID'] == cluster_id]['Epitope_ID'].tolist()
        cluster_to_epitopes[cluster_id] = epitope_ids
    
    print(f"\nLoaded cluster mapping: {len(cluster_to_epitopes)} clusters")
    
    # 2. Load protein mappings and reparse species
    mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/TOP50_protein_mapping.txt'
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    print(f"Loaded protein mappings: {len(mapping_df)} records")
    print("\nReparsing species with CORRECTED function...")
    
    mapping_df['Species_NEW'] = mapping_df['Full_Protein_Header'].apply(parse_species_from_protein_id)
    
    # 3. Calculate NEW species counts per cluster
    new_species_counts = []
    
    for cluster_id, epitope_ids in cluster_to_epitopes.items():
        cluster_mappings = mapping_df[mapping_df['Epitope_ID'].isin(epitope_ids)]
        
        species_set = set()
        for species in cluster_mappings['Species_NEW'].dropna():
            species_set.add(species)
        
        new_species_counts.append({
            'Cluster_ID': cluster_id,
            'Species_List_NEW': sorted(list(species_set)),
            'Num_Species_NEW': len(species_set),
            'Species_Names_NEW': '; '.join(sorted(list(species_set))) if species_set else 'None'
        })
    
    new_counts_df = pd.DataFrame(new_species_counts)
    
    # 4. Load OLD TOP50 metadata
    old_top50 = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt', sep='\t')
    
    print(f"Loaded OLD TOP50: {len(old_top50)} clusters")
    
    # 5. Merge and compare
    comparison = old_top50[['Cluster_ID', 'Rank', 'Representative_Sequence', 'Priority_Score', 
                            'Num_Species', 'Species_List', 'Methods_Used']].merge(
        new_counts_df,
        on='Cluster_ID',
        how='left'
    )
    
    comparison['Species_Count_Changed'] = comparison['Num_Species'] != comparison['Num_Species_NEW']
    comparison['Species_Diff'] = comparison['Num_Species_NEW'] - comparison['Num_Species']
    
    # 6. Summary statistics
    print("\n" + "="*80)
    print("IMPACT SUMMARY")
    print("="*80)
    
    print(f"\nTotal TOP50 clusters: {len(comparison)}")
    print(f"Clusters with CHANGED species counts: {comparison['Species_Count_Changed'].sum()}")
    print(f"Clusters with UNCHANGED species counts: {(~comparison['Species_Count_Changed']).sum()}")
    
    print(f"\nSpecies count changes:")
    diff_counts = comparison['Species_Diff'].value_counts().sort_index()
    for diff, count in diff_counts.items():
        if diff > 0:
            print(f"  +{diff} species: {count} clusters")
        elif diff < 0:
            print(f"  {diff} species: {count} clusters")
        else:
            print(f"  No change: {count} clusters")
    
    # 7. Show clusters with biggest changes
    print("\n" + "="*80)
    print("CLUSTERS WITH INCREASED SPECIES COUNTS (sorted by difference)")
    print("="*80)
    
    increased = comparison[comparison['Species_Diff'] > 0].sort_values('Species_Diff', ascending=False)
    
    if len(increased) > 0:
        print(f"\n{len(increased)} clusters gained species:\n")
        for _, row in increased.iterrows():
            print(f"Cluster {row['Cluster_ID']} (Rank {row['Rank']}): {row['Representative_Sequence']}")
            print(f"  OLD: {row['Num_Species']} species - NEW: {row['Num_Species_NEW']} species (+{row['Species_Diff']})")
            print(f"  Method: {row['Methods_Used']}")
            print(f"  OLD species: {row['Species_List']}")
            print(f"  NEW species: {row['Species_Names_NEW']}")
            print()
    else:
        print("No clusters gained species")
    
    # 8. Show clusters with decreased counts
    print("="*80)
    print("CLUSTERS WITH DECREASED SPECIES COUNTS")
    print("="*80)
    
    decreased = comparison[comparison['Species_Diff'] < 0].sort_values('Species_Diff')
    
    if len(decreased) > 0:
        print(f"\n{len(decreased)} clusters lost species:\n")
        for _, row in decreased.iterrows():
            print(f"Cluster {row['Cluster_ID']} (Rank {row['Rank']}): {row['Representative_Sequence']}")
            print(f"  OLD: {row['Num_Species']} species - NEW: {row['Num_Species_NEW']} species ({row['Species_Diff']})")
            print(f"  Method: {row['Methods_Used']}")
            print()
    else:
        print("No clusters lost species")
    
    # 9. Check if any non-TOP50 epitopes now have MORE species than current TOP50
    print("\n" + "="*80)
    print("CHECKING FOR EXCLUDED EPITOPES WITH HIGH SPECIES COUNTS")
    print("="*80)
    
    # Load all prioritized epitopes
    all_prioritized = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/PRIORITIZED_all_epitopes.txt', sep='\t')
    
    print(f"\nLoaded {len(all_prioritized)} total prioritized epitopes")
    
    # Get epitopes NOT in TOP50
    top50_epitopes = set(old_top50['Cluster_ID'])
    excluded_epitopes = all_prioritized[~all_prioritized['Epitope_ID'].isin(top50_epitopes)]
    
    print(f"Epitopes NOT in TOP50: {len(excluded_epitopes)}")
    
    # Current minimum species count in TOP50
    min_species_in_top50_old = old_top50['Num_Species'].min()
    min_species_in_top50_new = comparison['Num_Species_NEW'].min()
    
    print(f"\nMinimum species count in TOP50:")
    print(f"  OLD: {min_species_in_top50_old} species")
    print(f"  NEW: {min_species_in_top50_new} species")
    
    # Find excluded epitopes with high species counts (using NEW parsing)
    # We need to recalculate for ALL epitopes, not just TOP50
    print("\nRecalculating species counts for ALL epitopes with NEW parsing...")
    print("(This may take a moment...)")
    
    # Load all protein mappings (not just TOP50)
    all_mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/TOP50_protein_mapping.txt'
    # Note: This only has TOP50 mappings. We'd need the full mapping for all epitopes.
    # For now, let's just note this limitation
    
    print("\nNOTE: Full comparison requires protein mappings for ALL epitopes,")
    print("not just TOP50. This will be done in the full pipeline rebuild.")
    
    # 10. Save detailed comparison
    comparison_file = f"{output_dir}/TOP50_species_comparison_DETAILED.txt"
    comparison.to_csv(comparison_file, sep='\t', index=False)
    print(f"\n✓ Saved detailed comparison: {comparison_file}")
    
    # 11. Final recommendation
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    
    if comparison['Species_Count_Changed'].sum() > 0:
        print("\n⚠️  Species counts HAVE CHANGED for some clusters!")
        print("\nIMPACT:")
        print("  - Priority scores will change (they include species count)")
        print("  - Rankings may shift")
        print("  - Different epitopes might enter TOP50")
        print("\nRECOMMENDATION:")
        print("  → REBUILD the pipeline from prioritization onwards")
        print("  → This ensures consistency for publication")
    else:
        print("\n✓ No species count changes detected")
        print("  → Current TOP50 is still valid")
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()