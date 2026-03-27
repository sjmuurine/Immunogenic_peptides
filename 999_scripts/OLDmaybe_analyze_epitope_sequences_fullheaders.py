import pandas as pd

def analyze_epitope_sequences(input_file, output_file):
    """
    Groups BLAST results by epitope sequence similarity.
    Now works with full protein headers.
    
    Parameters:
    - input_file: epitope_blast_COMBINED_fullheaders_high_quality.txt
    - output_file: output with sequence clustering
    """
    
    # Read the BLAST results
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total high-quality hits: {len(df)}")
    print(f"Unique epitope IDs: {df['Epitope_ID'].nunique()}")
    print(f"Unique epitope sequences: {df['Epitope_Sequence'].nunique()}")
    
    # Group by epitope sequence (exact matches)
    sequence_groups = df.groupby('Epitope_Sequence')
    
    results = []
    
    for seq, group in sequence_groups:
        # Get info about this sequence group
        epitope_ids = group['Epitope_ID'].unique()
        num_ids = len(epitope_ids)
        num_hits_anatomia = len(group[group['Data_Source'] == 'ANATOMIA'])
        num_hits_literature = len(group[group['Data_Source'] == 'LITERATURE'])
        total_hits = len(group)
        
        # Get representative protein matches with FULL HEADERS
        top_anatomia = group[group['Data_Source'] == 'ANATOMIA'].nlargest(1, 'Percent_Identity')
        top_literature = group[group['Data_Source'] == 'LITERATURE'].nlargest(1, 'Percent_Identity')
        
        anatomia_protein = top_anatomia['Protein_Full_Header'].values[0] if len(top_anatomia) > 0 else "None"
        anatomia_identity = top_anatomia['Percent_Identity'].values[0] if len(top_anatomia) > 0 else 0
        
        literature_protein = top_literature['Protein_Full_Header'].values[0] if len(top_literature) > 0 else "None"
        literature_identity = top_literature['Percent_Identity'].values[0] if len(top_literature) > 0 else 0
        
        avg_identity = group['Percent_Identity'].mean()
        max_identity = group['Percent_Identity'].max()
        
        results.append({
            'Epitope_Sequence': seq,
            'Sequence_Length': len(seq),
            'Num_IEDB_IDs': num_ids,
            'IEDB_IDs': '; '.join(epitope_ids),
            'Total_Hits': total_hits,
            'ANATOMIA_Hits': num_hits_anatomia,
            'LITERATURE_Hits': num_hits_literature,
            'Avg_Identity': round(avg_identity, 2),
            'Max_Identity': round(max_identity, 2),
            'Best_ANATOMIA_Match': anatomia_protein,
            'Best_ANATOMIA_Identity': round(anatomia_identity, 2),
            'Best_LITERATURE_Match': literature_protein,
            'Best_LITERATURE_Identity': round(literature_identity, 2),
            'Found_In': 'Both' if (num_hits_anatomia > 0 and num_hits_literature > 0) 
                        else ('ANATOMIA' if num_hits_anatomia > 0 else 'LITERATURE')
        })
    
    # Create summary dataframe
    summary_df = pd.DataFrame(results)
    
    # Sort by total hits (most hits first)
    summary_df = summary_df.sort_values('Total_Hits', ascending=False)
    
    # Save to file
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print(f"{'='*60}")
    print(f"Unique epitope sequences: {len(summary_df)}")
    print(f"Sequences found in both datasets: {len(summary_df[summary_df['Found_In'] == 'Both'])}")
    print(f"Sequences only in ANATOMIA: {len(summary_df[summary_df['Found_In'] == 'ANATOMIA'])}")
    print(f"Sequences only in LITERATURE: {len(summary_df[summary_df['Found_In'] == 'LITERATURE'])}")
    print(f"\nTop 5 epitope sequences by total hits:")
    print(summary_df[['Epitope_Sequence', 'Total_Hits', 'ANATOMIA_Hits', 'LITERATURE_Hits', 'Found_In']].head())
    print(f"\n{'='*60}")
    print(f"Results saved to: {output_file}")
    print(f"{'='*60}")

# ============================================
# Run on Puhti or locally
# ============================================

input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
output_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_sequences_summary_fullheaders.txt'

analyze_epitope_sequences(input_file, output_file)
