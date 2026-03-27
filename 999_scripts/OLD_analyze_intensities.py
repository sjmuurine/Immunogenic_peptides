import pandas as pd
import os
import numpy as np

def analyze_epitope_intensities(epitope_files, intensity_file, output_dir):
    """
    Map epitopes to ANATOMIA protein intensities.
    Uses Mean relative intensity (protein-level) and OG median intensities.
    """
    
    print("="*70)
    print("EPITOPE INTENSITY ANALYSIS")
    print("="*70)
    
    # Load intensity data
    print(f"\nLoading intensity data: {intensity_file}")
    intensity_df = pd.read_csv(intensity_file, sep='\t')
    
    # Create protein ID -> intensity mapping
    intensity_map = {}
    for _, row in intensity_df.iterrows():
        protein_id = row['Protein_ID']
        intensity_map[protein_id] = {
            'species': row.get('Species', 'Unknown'),
            'orthogroup': row.get('Orthogroup', 'Unknown'),
            'protein_intensity': row.get('Protein_Mean_Relative_Intensity', np.nan),
            'og_intensity_all': row.get('OG_Median_Intensity_All', np.nan),
            'og_intensity_grampos': row.get('OG_Median_Intensity_GramPos', np.nan),
            'og_intensity_gramneg': row.get('OG_Median_Intensity_GramNeg', np.nan),
            'localization': row.get('Localization', 'Unknown')
        }
    
    print(f"Loaded intensity data for {len(intensity_map)} proteins")
    
    # Load protein mappings
    mapping_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    protein_mappings = {}
    
    for tier in ['TOP50', 'TIER1', 'TIER2']:
        mapping_file = f"{mapping_dir}/{tier}_protein_mapping.txt"
        try:
            protein_mappings[tier] = pd.read_csv(mapping_file, sep='\t')
            print(f"Loaded {tier} protein mappings: {len(protein_mappings[tier])} entries")
        except:
            protein_mappings[tier] = pd.DataFrame()
    
    # Process each epitope file
    for file_name, file_path in epitope_files.items():
        print(f"\n{'='*70}")
        print(f"Processing: {file_name}")
        print(f"{'='*70}")
        
        try:
            epitope_df = pd.read_csv(file_path, sep='\t')
        except:
            print(f"Could not load {file_path}")
            continue
        
        print(f"Loaded {len(epitope_df)} epitopes")
        
        # Get detailed cluster info
        detailed_file = f"/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/{file_name.replace('_FINAL', '')}_cluster_detailed_v2.txt"
        try:
            detailed_df = pd.read_csv(detailed_file, sep='\t')
        except:
            detailed_df = pd.DataFrame()
        
        results = []
        
        for _, epitope_row in epitope_df.iterrows():
            cluster_id = epitope_row['Cluster_ID']
            epitope_seq = epitope_row['Representative_Sequence']
            
            # Get all epitope IDs in this cluster
            if len(detailed_df) > 0:
                cluster_epitopes = detailed_df[detailed_df['Cluster_ID'] == cluster_id]
                epitope_ids = cluster_epitopes['Epitope_ID'].tolist()
            else:
                epitope_ids = [epitope_row.get('Representative_ID', cluster_id)]
            
            # Find ANATOMIA proteins for this epitope
            anatomia_proteins = []
            
            for epitope_id in epitope_ids:
                for tier, mapping_df in protein_mappings.items():
                    if len(mapping_df) > 0:
                        # Get ANATOMIA proteins only
                        matches = mapping_df[
                            (mapping_df['Epitope_ID'] == epitope_id) &
                            (mapping_df['Data_Source'] == 'ANATOMIA')
                        ]
                        
                        for _, match in matches.iterrows():
                            protein_id = match['Protein_ID']
                            
                            # Check if we have intensity data
                            if protein_id in intensity_map:
                                anatomia_proteins.append({
                                    'protein_id': protein_id,
                                    **intensity_map[protein_id]
                                })
            
            # Aggregate intensity data
            if anatomia_proteins:
                # Protein-level intensities
                protein_intensities = [p['protein_intensity'] for p in anatomia_proteins if pd.notna(p['protein_intensity'])]
                
                if protein_intensities:
                    sum_protein_intensity = sum(protein_intensities)
                    mean_protein_intensity = np.mean(protein_intensities)
                    max_protein_intensity = max(protein_intensities)
                    median_protein_intensity = np.median(protein_intensities)
                else:
                    sum_protein_intensity = 0
                    mean_protein_intensity = 0
                    max_protein_intensity = 0
                    median_protein_intensity = 0
                
                # OG-level intensities
                og_intensities_all = [p['og_intensity_all'] for p in anatomia_proteins if pd.notna(p['og_intensity_all'])]
                
                if og_intensities_all:
                    mean_og_intensity = np.mean(og_intensities_all)
                    max_og_intensity = max(og_intensities_all)
                else:
                    mean_og_intensity = 0
                    max_og_intensity = 0
                
                # Get unique species and OGs with intensity
                species_with_intensity = set(p['species'] for p in anatomia_proteins if p['species'] != 'Unknown')
                ogs_with_intensity = set(p['orthogroup'] for p in anatomia_proteins if p['orthogroup'] != 'Unknown')
                
                protein_ids = [p['protein_id'] for p in anatomia_proteins[:10]]
                
                num_with_protein_intensity = len(protein_intensities)
                num_with_og_intensity = len(og_intensities_all)
            else:
                sum_protein_intensity = 0
                mean_protein_intensity = 0
                max_protein_intensity = 0
                median_protein_intensity = 0
                mean_og_intensity = 0
                max_og_intensity = 0
                species_with_intensity = set()
                ogs_with_intensity = set()
                protein_ids = []
                num_with_protein_intensity = 0
                num_with_og_intensity = 0
            
            results.append({
                'Cluster_ID': cluster_id,
                'Rank': epitope_row.get('Rank', ''),
                'Epitope_Sequence': epitope_seq,
                'Priority_Score': epitope_row.get('Priority_Score', ''),
                'Methods_Used': epitope_row.get('Methods_Used', ''),
                'Suitable_For_Fishing': epitope_row.get('Suitable_For_Fishing', ''),
                
                # Protein-level intensity data
                'Num_ANATOMIA_Proteins': len(anatomia_proteins),
                'Num_With_Protein_Intensity': num_with_protein_intensity,
                'Sum_Protein_Intensity': sum_protein_intensity,
                'Mean_Protein_Intensity': mean_protein_intensity,
                'Median_Protein_Intensity': median_protein_intensity,
                'Max_Protein_Intensity': max_protein_intensity,
                
                # OG-level intensity data
                'Num_With_OG_Intensity': num_with_og_intensity,
                'Mean_OG_Intensity': mean_og_intensity,
                'Max_OG_Intensity': max_og_intensity,
                
                # Species and OG info
                'Species_With_Intensity': '; '.join(sorted(species_with_intensity)) if species_with_intensity else 'None',
                'Num_Species_With_Intensity': len(species_with_intensity),
                'OGs_With_Intensity': '; '.join(sorted(ogs_with_intensity)) if ogs_with_intensity else 'None',
                'Num_OGs_With_Intensity': len(ogs_with_intensity),
                
                'Example_Proteins': '; '.join(protein_ids[:5]) if protein_ids else 'None',
                
                # Original metadata
                'Total_Predicted_Species': epitope_row.get('Num_Species', ''),
                'Total_ANATOMIA_Hits': epitope_row.get('ANATOMIA_Hits', ''),
            })
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Sort by mean protein intensity (most relevant metric)
        results_df = results_df.sort_values('Mean_Protein_Intensity', ascending=False)
        
        # Save
        output_file = f"{output_dir}/{file_name}_with_intensities.txt"
        results_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\n✓ Saved: {output_file}")
        
        # Statistics
        with_protein_intensity = len(results_df[results_df['Num_With_Protein_Intensity'] > 0])
        with_og_intensity = len(results_df[results_df['Num_With_OG_Intensity'] > 0])
        with_any = len(results_df[(results_df['Num_With_Protein_Intensity'] > 0) | (results_df['Num_With_OG_Intensity'] > 0)])
        
        print(f"\nResults:")
        print(f"  Epitopes with protein-level intensity: {with_protein_intensity}/{len(results_df)}")
        print(f"  Epitopes with OG-level intensity: {with_og_intensity}/{len(results_df)}")
        print(f"  Epitopes with any intensity data: {with_any}/{len(results_df)}")
        
        if with_protein_intensity > 0:
            print(f"\nTop 5 by mean protein intensity:")
            top5 = results_df[results_df['Num_With_Protein_Intensity'] > 0].head(5)
            for _, row in top5.iterrows():
                print(f"  {row['Epitope_Sequence']:20} - {row['Num_With_Protein_Intensity']} proteins, mean={row['Mean_Protein_Intensity']:.4f}")
    
    print("\n" + "="*70)
    print("INTENSITY ANALYSIS COMPLETE")
    print("="*70)

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_intensity_analysis'
os.makedirs(output_dir, exist_ok=True)

epitope_files = {
    'TOP50_FINAL': '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt',
    'TIER1_FINAL': '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TIER1_FINAL_with_metadata.txt',
    'TIER2_FINAL': '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TIER2_FINAL_with_metadata.txt',
    'COMBINED_TOP_TIERS_FINAL': '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/COMBINED_TOP_TIERS_FINAL_with_metadata.txt',
    'COMBINED_FISHABLE': '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/COMBINED_TOP_TIERS_FISHABLE_epitopes.txt'
}

intensity_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'

analyze_epitope_intensities(epitope_files, intensity_file, output_dir)