import pandas as pd

print("="*70)
print("CHECKING FILE STRUCTURES FOR METHOD COMPARISON")
print("="*70)

# 1. IEDB BLAST results
print("\n1. IEDB BLAST FILE:")
iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
iedb_df = pd.read_csv(iedb_file, sep='\t', nrows=5)
print(f"   Rows (sample): {len(iedb_df)}")
print(f"   Columns: {iedb_df.columns.tolist()}")
print(f"\n   Sample data:")
print(iedb_df[['Epitope_Sequence', 'Protein_Full_Header', 'Data_Source']].head(3))

# Count total
iedb_full = pd.read_csv(iedb_file, sep='\t')
print(f"\n   Total IEDB records: {len(iedb_full)}")
print(f"   Unique epitope sequences: {iedb_full['Epitope_Sequence'].nunique()}")

# 2. K-mer EXACT results
print("\n\n2. K-MER EXACT FILE:")
exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
exact_df = pd.read_csv(exact_file, sep='\t', nrows=5)
print(f"   Rows (sample): {len(exact_df)}")
print(f"   Columns: {exact_df.columns.tolist()}")
print(f"\n   Sample data:")
print(exact_df.head(3))

exact_full = pd.read_csv(exact_file, sep='\t')
print(f"\n   Total exact motifs: {len(exact_full)}")

# 3. K-mer FUZZY results
print("\n\n3. K-MER FUZZY FILE:")
fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
fuzzy_df = pd.read_csv(fuzzy_file, nrows=5)
print(f"   Rows (sample): {len(fuzzy_df)}")
print(f"   Columns: {fuzzy_df.columns.tolist()}")
print(f"\n   Sample data:")
print(fuzzy_df.head(3))

fuzzy_full = pd.read_csv(fuzzy_file)
print(f"\n   Total fuzzy records: {len(fuzzy_full)}")
print(f"   Unique consensus sequences: {fuzzy_full['Consensus'].nunique()}")

print("\n" + "="*70)
print("QUESTIONS:")
print("="*70)
print("\n1. For IEDB: Use 'Epitope_Sequence' column?")
print("2. For EXACT: Use 'Motif' column?")
print("3. For FUZZY: Use 'Consensus' column?")
print("4. Should we match EXACT sequences (character-by-character)?")
print("5. Or should FUZZY matches count too (with similarity threshold)?")
print("\n" + "="*70)