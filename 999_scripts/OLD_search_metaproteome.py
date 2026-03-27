import pandas as pd
import os

def search_epitopes_in_metaproteome(epitope_files, metaproteome_fasta, output_dir):
    """
    Search epitopes in metaproteome proteins using exact string matching.
    """
    
    print("="*70)
    print("SEARCHING EPITOPES IN METAPROTEOME")
    print("="*70)
    
    # Load metaproteome sequences
    print(f"\nLoading metaproteome: {metaproteome_fasta}")
    metaproteome_proteins = {}
    
    current_id = None
    current_seq = []
    
    with open(metaproteome_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous protein
                if current_id and current_seq:
                    metaproteome_proteins[current_id]['sequence'] = ''.join(current_seq)
                
                # Parse new header: >ProteinID|Intensity|Genus|Gram|Name
                header = line[1:]
                parts = header.split('|')
                
                protein_id = parts[0] if len(parts) > 0 else 'Unknown'
                
                # Initialize with defaults
                metaproteome_proteins[protein_id] = {
                    'intensity': 0.0,
                    'genus': 'Unknown',
                    'gram': 'Unknown',
                    'name': 'Unknown',
                    'sequence': ''
                }
                
                # Parse available fields
                try:
                    if len(parts) > 1:
                        metaproteome_proteins[protein_id]['intensity'] = float(parts[1])
                    if len(parts) > 2:
                        metaproteome_proteins[protein_id]['genus'] = parts[2]
                    if len(parts) > 3:
                        metaproteome_proteins[protein_id]['gram'] = parts[3]
                    if len(parts) > 4:
                        metaproteome_proteins[protein_id]['name'] = parts[4]
                except:
                    pass  # Keep defaults
                
                current_id = protein_id
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last protein
        if current_id and current_seq:
            metaproteome_proteins[current_id]['sequence'] = ''.join(current_seq)
    
    print(f"Loaded {len(metaproteome_proteins)} metaproteome proteins")
    
    # Process each epitope file
    for epitope_file_name, epitope_file_path in epitope_files.items():
        print(f"\n{'='*70}")
        print(f"Processing: {epitope_file_name}")
        print(f"{'='*70}")
        
        try:
            epitope_df = pd.read_csv(epitope_file_path, sep='\t')
        except:
            print(f"Could not load {epitope_file_path}, skipping...")
            continue
        
        print(f"Loaded {len(epitope_df)} epitope clusters")
        
        results = []
        
        for _, epitope_row in epitope_df.iterrows():
            epitope_seq = epitope_row['Representative_Sequence']
            cluster_id = epitope_row['Cluster_ID']
            
            # Search for epitope in all metaproteome proteins
            found_in = []
            total_intensity = 0
            
            for protein_id, protein_data in metaproteome_proteins.items():
                protein_seq = protein_data.get('sequence', '')
                
                # Exact string match
                if epitope_seq in protein_seq:
                    found_in.append({
                        'protein_id': protein_id,
                        'intensity': protein_data.get('intensity', 0),
                        'genus': protein_data.get('genus', 'Unknown'),
                        'gram': protein_data.get('gram', 'Unknown'),
                        'name': protein_data.get('name', 'Unknown')
                    })
                    total_intensity += protein_data.get('intensity', 0)
            
            # Compile results
            results.append({
                'Cluster_ID': cluster_id,
                'Rank': epitope_row.get('Rank', ''),
                'Epitope_Sequence': epitope_seq,
                'Sequence_Length': len(epitope_seq),
                'Priority_Score': epitope_row.get('Priority_Score', ''),
                'Methods_Used': epitope_row.get('Methods_Used', ''),
                'IEDB_Links': epitope_row.get('IEDB_Links', ''),
                'Predicted_Species_List': epitope_row.get('Species_List', ''),
                'Predicted_Num_Species': epitope_row.get('Num_Species', ''),
                'Suitable_For_Fishing': epitope_row.get('Suitable_For_Fishing', ''),
                
                # Metaproteome findings
                'Found_In_Metaproteome': 'Yes' if len(found_in) > 0 else 'No',
                'Num_Metaproteome_Proteins': len(found_in),
                'Total_Intensity': total_intensity,
                'Mean_Intensity': total_intensity / len(found_in) if len(found_in) > 0 else 0,
                
                # Protein details (show up to 10)
                'Metaproteome_Protein_IDs': '; '.join([p['protein_id'] for p in found_in[:10]]) if found_in else 'None',
                'Metaproteome_Genera': '; '.join(sorted(set([p['genus'] for p in found_in]))) if found_in else 'None',
                'Metaproteome_Gram': '; '.join(sorted(set([p['gram'] for p in found_in]))) if found_in else 'None',
            })
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Sort by whether found, then by intensity
        results_df['Found_Sort'] = results_df['Found_In_Metaproteome'].map({'Yes': 0, 'No': 1})
        results_df = results_df.sort_values(['Found_Sort', 'Total_Intensity'], ascending=[True, False])
        results_df = results_df.drop('Found_Sort', axis=1)
        
        # Save
        output_file = f"{output_dir}/{epitope_file_name}_metaproteome_search.txt"
        results_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"\n✓ Saved: {output_file}")
        
        # Statistics
        found = len(results_df[results_df['Found_In_Metaproteome'] == 'Yes'])
        not_found = len(results_df[results_df['Found_In_Metaproteome'] == 'No'])
        
        print(f"\nResults:")
        print(f"  Found in metaproteome: {found}/{len(results_df)} ({100*found/len(results_df):.1f}%)")
        print(f"  Not found: {not_found}/{len(results_df)}")
        
        if found > 0:
            top_found = results_df[results_df['Found_In_Metaproteome'] == 'Yes'].head(5)
            print(f"\nTop 5 found epitopes by intensity:")
            for _, row in top_found.iterrows():
                print(f"  {row['Epitope_Sequence']} - {row['Num_Metaproteome_Proteins']} proteins, intensity={row['Total_Intensity']:.2e}")
    
    print("\n" + "="*70)
    print("SEARCH COMPLETE")
    print("="*70)

# Run
cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
metaproteome_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/metaproteome_sequences.fa'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Define epitope files to search
epitope_files = {
    'TOP50_FINAL': f"{cluster_dir}/TOP50_FINAL_with_metadata.txt",
    'TIER1_FINAL': f"{cluster_dir}/TIER1_FINAL_with_metadata.txt",
    'TIER2_FINAL': f"{cluster_dir}/TIER2_FINAL_with_metadata.txt",
    'COMBINED_TOP_TIERS_FINAL': f"{cluster_dir}/COMBINED_TOP_TIERS_FINAL_with_metadata.txt",
    'TOP50_FISHABLE': f"{cluster_dir}/TOP50_FISHABLE_epitopes.txt",
    'COMBINED_FISHABLE': f"{cluster_dir}/COMBINED_TOP_TIERS_FISHABLE_epitopes.txt"
}

search_epitopes_in_metaproteome(epitope_files, metaproteome_fasta, output_dir)