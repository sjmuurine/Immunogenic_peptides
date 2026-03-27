import pandas as pd

def parse_species_from_protein_id(protein_id):
    """Parse species from protein ID - CORRECTED VERSION."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
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
    
    if 'LITERATURE' in protein_id:
        for pattern, species in species_patterns.items():
            if pattern.lower() in protein_id.lower():
                return species
        return None
    
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

print("="*80)
print("INVESTIGATING FUZZY EPITOPES WITH 16 SPECIES")
print("="*80)

# Load fuzzy species distribution
fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/fuzzy_species_distribution.txt'
fuzzy_df = pd.read_csv(fuzzy_file, sep='\t')

print(f"\nTotal fuzzy epitopes: {len(fuzzy_df)}")

# Find epitopes with 16 species
high_species = fuzzy_df[fuzzy_df['Num_Species'] >= 16].sort_values('Num_Species', ascending=False)

print(f"\nFuzzy epitopes with 16+ species: {len(high_species)}")
print("\nThese epitopes:")
for _, row in high_species.iterrows():
    print(f"\n  Sequence: {row['Epitope_Sequence']}")
    print(f"  Species count: {row['Num_Species']}")
    print(f"  Species: {row['Species_Names']}")

# Now check: are these in the BACTERIA_SPECIFIC file (passed human filter)?
print("\n" + "="*80)
print("CHECKING IF THESE EPITOPES PASSED HUMAN FILTER")
print("="*80)

bacteria_specific = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/BACTERIA_SPECIFIC_epitopes.txt', sep='\t')

print(f"\nTotal bacteria-specific epitopes: {len(bacteria_specific)}")
print(f"Columns in bacteria_specific file: {bacteria_specific.columns.tolist()}")

# Find the sequence column
seq_col = None
for col in bacteria_specific.columns:
    if 'sequence' in col.lower() or 'epitope' in col.lower():
        seq_col = col
        break

if seq_col:
    print(f"Using sequence column: '{seq_col}'")
    
    for _, row in high_species.iterrows():
        seq = row['Epitope_Sequence']
        in_bacteria = bacteria_specific[bacteria_specific[seq_col] == seq]
        
        if len(in_bacteria) > 0:
            print(f"\n✓ {seq}")
            print(f"  Status: PASSED human filter")
            print(f"  Found in bacteria_specific file")
            available_cols = in_bacteria.columns.tolist()
            if 'Method' in available_cols:
                print(f"  Method: {in_bacteria['Method'].values[0]}")
            if 'Max_Human_Identity' in available_cols:
                print(f"  Max human identity: {in_bacteria['Max_Human_Identity'].values[0]:.1f}%")
        else:
            print(f"\n✗ {seq}")
            print(f"  Status: FAILED human filter (excluded)")
else:
    print("ERROR: Could not find sequence column!")
    print("Please check the file manually")

# Check if they're in the PRIORITIZED file
print("\n" + "="*80)
print("CHECKING IF THESE EPITOPES WERE PRIORITIZED")
print("="*80)

prioritized = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/PRIORITIZED_all_epitopes.txt', sep='\t')

print(f"\nTotal prioritized epitopes: {len(prioritized)}")
print(f"Columns in prioritized file: {prioritized.columns.tolist()}")

# Find sequence column
seq_col_prior = None
for col in prioritized.columns:
    if 'sequence' in col.lower() or 'epitope' in col.lower():
        seq_col_prior = col
        break

if seq_col_prior:
    print(f"Using sequence column: '{seq_col_prior}'")
    
    for _, row in high_species.iterrows():
        seq = row['Epitope_Sequence']
        in_prioritized = prioritized[prioritized[seq_col_prior] == seq]
        
        if len(in_prioritized) > 0:
            print(f"\n✓ {seq}")
            print(f"  Status: WAS PRIORITIZED")
            if 'Priority_Score' in in_prioritized.columns:
                print(f"  Priority score: {in_prioritized['Priority_Score'].values[0]}")
            if 'Methods_Used' in in_prioritized.columns:
                print(f"  Methods: {in_prioritized['Methods_Used'].values[0]}")
        else:
            print(f"\n✗ {seq}")
            print(f"  Status: NOT in prioritized list")

# Check TOP50
print("\n" + "="*80)
print("CHECKING IF THESE EPITOPES MADE IT TO TOP50")
print("="*80)

top50 = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt', sep='\t')

for _, row in high_species.iterrows():
    seq = row['Epitope_Sequence']
    in_top50 = top50[top50['Representative_Sequence'] == seq]
    
    if len(in_top50) > 0:
        print(f"\n✓ {seq}")
        print(f"  Status: IN TOP50")
        print(f"  Rank: {in_top50['Rank'].values[0]}")
    else:
        print(f"\n✗ {seq}")
        print(f"  Status: NOT in TOP50")

# Compare with current TOP50 minimum
print("\n" + "="*80)
print("COMPARISON WITH CURRENT TOP50")
print("="*80)

print(f"\nHighest species count in TOP50: {top50['Num_Species'].max()}")
print(f"Lowest species count in TOP50: {top50['Num_Species'].min()}")
print(f"Median species count in TOP50: {top50['Num_Species'].median()}")

print(f"\nFuzzy epitopes with 16 species would outrank (by species count alone):")
print(f"  {(top50['Num_Species'] < 16).sum()} epitopes in current TOP50")

# Check priority scoring by method
print("\n" + "="*80)
print("PRIORITY SCORING INVESTIGATION")
print("="*80)

print("\nLet's check if fuzzy k-mer gets lower priority scores...")

if seq_col_prior and 'Methods_Used' in prioritized.columns and 'Priority_Score' in prioritized.columns:
    fuzzy_only = prioritized[prioritized['Methods_Used'].str.contains('FUZZY', na=False)]
    exact_only = prioritized[prioritized['Methods_Used'].str.contains('EXACT', na=False) & ~prioritized['Methods_Used'].str.contains('FUZZY', na=False)]
    
    print(f"\nFUZZY (any combination) in PRIORITIZED:")
    print(f"  Count: {len(fuzzy_only)}")
    if len(fuzzy_only) > 0:
        print(f"  Average priority score: {fuzzy_only['Priority_Score'].mean():.1f}")
        print(f"  Max priority score: {fuzzy_only['Priority_Score'].max():.1f}")
        print(f"  Min priority score: {fuzzy_only['Priority_Score'].min():.1f}")
    
    print(f"\nEXACT (without fuzzy) in PRIORITIZED:")
    print(f"  Count: {len(exact_only)}")
    if len(exact_only) > 0:
        print(f"  Average priority score: {exact_only['Priority_Score'].mean():.1f}")
        print(f"  Max priority score: {exact_only['Priority_Score'].max():.1f}")
        print(f"  Min priority score: {exact_only['Priority_Score'].min():.1f}")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)

print("\nPossible reasons fuzzy epitopes with 16 species didn't make TOP50:")
print("  1. Failed human similarity filter")
print("  2. Lower priority score due to method weighting")
print("  3. Not considered bacteria-specific enough")
print("  4. Lower total hits count")
print("\nCheck the output above to see which one applies!")

print("\n" + "="*80)