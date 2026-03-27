import pandas as pd
import os

def analyze_fuzzy_species_distribution(output_dir):
    """
    Analyze species distribution of fuzzy k-mer epitopes.
    """
    
    print("="*70)
    print("ANALYZING FUZZY K-MER SPECIES DISTRIBUTION")
    print("="*70)
    
    # Load fuzzy epitope data
    fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
    
    print(f"\nLoading: {fuzzy_file}")
    fuzzy_df = pd.read_csv(fuzzy_file)
    
    print(f"Total fuzzy epitope records: {len(fuzzy_df)}")
    print(f"Unique consensus sequences: {fuzzy_df['Consensus'].nunique()}")
    
    # Load ANATOMIA protein info
    anatomia_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/ANATOMIA_bEVs_all_sequences.fa'
    
    print(f"\nLoading ANATOMIA protein headers from: {anatomia_file}")
    
    # Parse FASTA headers
    protein_headers = {}
    with open(anatomia_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Format: >WP_012420089.1 Species_OG_Conservation_Localization
                header = line.strip()[1:]  # Remove '>'
                parts = header.split(' ', 1)
                if len(parts) == 2:
                    protein_id = parts[0]
                    full_header = header
                    protein_headers[protein_id] = full_header
    
    print(f"Loaded {len(protein_headers)} protein headers")
    
    # Map Source_Protein to full headers
    fuzzy_df['Full_Header'] = fuzzy_df['Source_Protein'].map(protein_headers)
    
    print(f"Mapped headers: {fuzzy_df['Full_Header'].notna().sum()} / {len(fuzzy_df)}")
    
    # Parse species from headers
    def parse_species_from_header(header):
        if pd.isna(header):
            return None
        
        header_str = str(header)
        
        # ANATOMIA format: "ProteinID Species_OG_Conservation_Localization"
        if ' ' not in header_str:
            return None
        
        parts = header_str.split(' ', 1)[1].split('_')
        
        # Find OG position
        og_idx = None
        for i, part in enumerate(parts):
            if part.startswith('OG') or part == 'unique':
                og_idx = i
                break
        
        if og_idx and og_idx > 0:
            return '_'.join(parts[:og_idx])
        
        return None
    
    # Parse species
    print("\nParsing species from headers...")
    fuzzy_df['Species'] = fuzzy_df['Full_Header'].apply(parse_species_from_header)
    
    print(f"Species parsed: {fuzzy_df['Species'].notna().sum()}")
    print(f"Unique species: {fuzzy_df['Species'].nunique()}")
    
    # Filter to records with species
    fuzzy_with_species = fuzzy_df[fuzzy_df['Species'].notna()].copy()
    
    print(f"\nRecords with species info: {len(fuzzy_with_species)}")
    
    # Group by epitope sequence and count species
    print("\nAnalyzing species distribution per epitope...")
    
    sequence_species = fuzzy_with_species.groupby('Consensus').agg({
        'Species': lambda x: sorted(list(x.unique())),
        'Source_Protein': 'count'
    }).reset_index()
    
    sequence_species.columns = ['Epitope_Sequence', 'Species_List', 'Total_Hits']
    sequence_species['Num_Species'] = sequence_species['Species_List'].apply(len)
    sequence_species['Species_Names'] = sequence_species['Species_List'].apply(lambda x: '; '.join(x))
    
    # Sort by number of species (ascending to see the narrow ones first)
    sequence_species = sequence_species.sort_values('Num_Species', ascending=True)
    
    # Save summary
    output_file = f"{output_dir}/fuzzy_species_distribution.txt"
    sequence_species[['Epitope_Sequence', 'Num_Species', 'Total_Hits', 'Species_Names']].to_csv(
        output_file, sep='\t', index=False
    )
    
    print(f"\n✓ Saved: {output_file}")
    
    # Statistics
    print("\n" + "="*70)
    print("SPECIES DISTRIBUTION STATISTICS")
    print("="*70)
    
    print(f"\nTotal unique fuzzy epitopes: {len(sequence_species)}")
    print(f"\nSpecies count distribution:")
    species_counts = sequence_species['Num_Species'].value_counts().sort_index()
    for num_species, count in species_counts.items():
        print(f"  {num_species} species: {count} epitopes")
    
    print(f"\nEpitopes found in only 1 species: {(sequence_species['Num_Species'] == 1).sum()}")
    print(f"Epitopes found in 2-3 species: {((sequence_species['Num_Species'] >= 2) & (sequence_species['Num_Species'] <= 3)).sum()}")
    print(f"Epitopes found in 4+ species: {(sequence_species['Num_Species'] >= 4).sum()}")
    
    if len(sequence_species) > 0:
        print(f"\nTop 10 epitopes by species count:")
        top10 = sequence_species.nlargest(10, 'Num_Species')
        for _, row in top10.iterrows():
            print(f"  {row['Epitope_Sequence']:22} - {row['Num_Species']:2} species: {row['Species_Names']}")
        
        print(f"\nBottom 10 epitopes (narrowest distribution):")
        bottom10 = sequence_species.nsmallest(10, 'Num_Species')
        for _, row in bottom10.iterrows():
            print(f"  {row['Epitope_Sequence']:22} - {row['Num_Species']:2} species: {row['Species_Names']}")
    
    # Create detailed per-sequence breakdown
    print("\nCreating detailed breakdown...")
    
    detailed_list = []
    for seq in sequence_species['Epitope_Sequence'].unique():
        seq_data = fuzzy_with_species[fuzzy_with_species['Consensus'] == seq]
        
        for species in seq_data['Species'].unique():
            species_proteins = seq_data[seq_data['Species'] == species]
            
            detailed_list.append({
                'Epitope_Sequence': seq,
                'Species': species,
                'Num_Proteins': len(species_proteins),
                'Num_Unique_Proteins': species_proteins['Source_Protein'].nunique(),
                'Example_Proteins': '; '.join(species_proteins['Source_Protein'].unique()[:3].tolist())
            })
    
    detailed_df = pd.DataFrame(detailed_list)
    detailed_df = detailed_df.sort_values(['Epitope_Sequence', 'Species'])
    
    detailed_file = f"{output_dir}/fuzzy_species_detailed.txt"
    detailed_df.to_csv(detailed_file, sep='\t', index=False)
    
    print(f"✓ Saved detailed breakdown: {detailed_file}")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
os.makedirs(output_dir, exist_ok=True)

analyze_fuzzy_species_distribution(output_dir)