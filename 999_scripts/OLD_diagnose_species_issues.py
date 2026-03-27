import pandas as pd

print("="*70)
print("DIAGNOSTIC: Finding problematic species parsing")
print("="*70)

# Check EXACT k-mer for LITERATURE_Unknown
print("\n1. CHECKING EXACT K-MER LITERATURE_Unknown:")
exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
exact_df = pd.read_csv(exact_file, sep='\t', nrows=100)

for idx, row in exact_df.iterrows():
    seq_ids = str(row['Sequence_IDs'])
    protein_ids = [pid.strip() for pid in seq_ids.split(',')]
    
    for protein_id in protein_ids:
        if 'LITERATURE' in protein_id:
            # Try to find species
            found_species = False
            species_patterns = [
                'fragilis', 'nucleatum', 'muciniphila', 'prausnitzii', 'reuteri',
                'cloacae', 'faecalis', 'salivarius', 'longum', 'plantarum',
                'casei', 'sanguinis', 'magna', 'russellii', 'freudenreichii', 'fermentum', 'adiacens'
            ]
            
            for pattern in species_patterns:
                if pattern in protein_id.lower():
                    found_species = True
                    break
            
            if not found_species:
                print(f"  UNMATCHED: {protein_id}")
                if idx > 5:  # Just show a few examples
                    break
    if idx > 5:
        break

# Check IEDB BLAST
print("\n2. CHECKING IEDB BLAST DATA_SOURCE:")
iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
iedb_df = pd.read_csv(iedb_file, sep='\t')

print(f"Total IEDB records: {len(iedb_df)}")
print(f"\nData source counts:")
print(iedb_df['Data_Source'].value_counts())

print(f"\nSample LITERATURE records:")
lit_samples = iedb_df[iedb_df['Data_Source'] == 'LITERATURE'].head(10)
print(lit_samples[['Protein_Full_Header', 'Data_Source']])

print("\n" + "="*70)