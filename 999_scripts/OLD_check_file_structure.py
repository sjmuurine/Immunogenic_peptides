import pandas as pd

# Check exact k-mer file
exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
exact_df = pd.read_csv(exact_file, sep='\t', nrows=5)
print("EXACT K-MER FILE:")
print(f"Columns: {exact_df.columns.tolist()}")
print(exact_df.head())

print("\n" + "="*70 + "\n")

# Check IEDB BLAST file
iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
iedb_df = pd.read_csv(iedb_file, sep='\t', nrows=5)
print("IEDB BLAST FILE:")
print(f"Columns: {iedb_df.columns.tolist()}")
print(iedb_df.head())