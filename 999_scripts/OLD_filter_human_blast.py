import pandas as pd

output_dir = "/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons"
temp_dir = f"{output_dir}/temp"

# Define column names
columns = ['Epitope_ID', 'Human_Protein', 'Percent_Identity', 'Alignment_Length', 
           'Epitope_Length', 'Human_Protein_Length', 'Query_Start', 'Query_End',
           'Subject_Start', 'Subject_End', 'E_value', 'Bit_Score', 'Query_Coverage',
           'Epitope_Sequence', 'Human_Sequence']

# Read BLAST results
print("Loading BLAST results...")
blast_df = pd.read_csv(f"{output_dir}/epitopes_vs_human_ALL_hits.txt", 
                       sep='\t', names=columns)

print(f"Total BLAST hits: {len(blast_df)}")
print(f"Unique epitopes with human hits: {blast_df['Epitope_ID'].nunique()}")

# Calculate alignment coverage
blast_df['Alignment_Coverage'] = (blast_df['Alignment_Length'] / blast_df['Epitope_Length']) * 100

# Filter criteria:
# Identity >= 75% AND alignment coverage >= 70% = TOO SIMILAR to human (reject)
high_human_similarity = blast_df[
    (blast_df['Percent_Identity'] >= 75) & 
    (blast_df['Alignment_Coverage'] >= 70)
]

# Get unique epitopes with high human similarity (to exclude)
epitopes_too_similar_to_human = set(high_human_similarity['Epitope_ID'].unique())

print(f"\nEpitopes with HIGH human similarity (≥75% identity, ≥70% coverage): {len(epitopes_too_similar_to_human)}")

# Load all epitopes
print("\nLoading all epitopes...")
all_epitopes = pd.read_csv(f"{temp_dir}/epitope_source_mapping.txt", sep='\t')

print(f"Total epitopes tested: {len(all_epitopes)}")

# Mark bacteria-specific epitopes
all_epitopes['Bacteria_Specific'] = ~all_epitopes['Epitope_ID'].isin(epitopes_too_similar_to_human)
bacteria_specific = all_epitopes[all_epitopes['Bacteria_Specific'] == True]

print(f"Bacteria-specific epitopes (LOW human similarity): {len(bacteria_specific)}")

# Add human similarity info
def get_best_human_match(epitope_id):
    matches = blast_df[blast_df['Epitope_ID'] == epitope_id]
    if len(matches) == 0:
        return 'No_hits', 0, 0
    best = matches.nlargest(1, 'Percent_Identity').iloc[0]
    return best['Human_Protein'][:100], best['Percent_Identity'], best['Alignment_Coverage']

print("\nAdding human similarity information...")
all_epitopes[['Best_Human_Match', 'Max_Human_Identity', 'Max_Human_Coverage']] = \
    all_epitopes['Epitope_ID'].apply(lambda x: pd.Series(get_best_human_match(x)))

# Save results
print("\nSaving results...")
all_epitopes.to_csv(f"{output_dir}/all_epitopes_with_human_similarity.txt", sep='\t', index=False)
bacteria_specific.to_csv(f"{output_dir}/BACTERIA_SPECIFIC_epitopes.txt", sep='\t', index=False)

# Save high human similarity epitopes (for reference)
epitopes_with_high_human_sim = all_epitopes[all_epitopes['Bacteria_Specific'] == False]
epitopes_with_high_human_sim.to_csv(f"{output_dir}/HIGH_HUMAN_SIMILARITY_epitopes_excluded.txt", sep='\t', index=False)

print("\n" + "="*60)
print("FILTERING COMPLETE")
print("="*60)
print(f"✓ All epitopes with human similarity info saved")
print(f"✓ Bacteria-specific epitopes: {len(bacteria_specific)}")
print(f"✓ High human similarity (excluded): {len(epitopes_with_high_human_sim)}")
print("="*60)