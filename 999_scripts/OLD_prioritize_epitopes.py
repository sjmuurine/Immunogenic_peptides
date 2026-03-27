import pandas as pd
import numpy as np
from collections import defaultdict

def get_sequence_metadata(seq, iedb_df, exact_df, fuzzy_detailed_df):
    """
    Get detailed metadata about where this sequence appears.
    """
    metadata = {
        'ANATOMIA_hits': 0,
        'LITERATURE_hits': 0,
        'Total_proteins': 0,
        'Orthogroups': set(),
        'Species': set(),
        'IEDB_IDs': set()
    }
    
    # Check IEDB data
    iedb_matches = iedb_df[iedb_df['Epitope_Sequence'] == seq]
    if len(iedb_matches) > 0:
        metadata['ANATOMIA_hits'] = len(iedb_matches[iedb_matches['Data_Source'] == 'ANATOMIA'])
        metadata['LITERATURE_hits'] = len(iedb_matches[iedb_matches['Data_Source'] == 'LITERATURE'])
        metadata['IEDB_IDs'] = set(iedb_matches['Epitope_ID'].unique())
    
    # Check k-mer exact data
    exact_matches = exact_df[exact_df['Motif'] == seq]
    if len(exact_matches) > 0:
        # Parse sequence IDs to get protein info
        seq_ids = str(exact_matches.iloc[0]['Sequence_IDs'])
        proteins = [p.strip() for p in seq_ids.split(',')]
        
        for protein in proteins:
            metadata['Total_proteins'] += 1
            # Extract species and OG from protein ID
            # Format: WP_xxxxx Species_OG_Conservation_Location
            parts = protein.split('_')
            if len(parts) >= 2:
                metadata['Species'].add(parts[1] if not parts[0].startswith('LITERATURE') else 'LITERATURE')
            if 'OG' in protein:
                og = [p for p in parts if p.startswith('OG')]
                if og:
                    metadata['Orthogroups'].add(og[0])
    
    # Check fuzzy data
    fuzzy_matches = fuzzy_detailed_df[fuzzy_detailed_df['Consensus'] == seq]
    if len(fuzzy_matches) > 0:
        for _, row in fuzzy_matches.iterrows():
            protein = str(row['Source_Protein'])
            if 'LITERATURE' in protein:
                metadata['LITERATURE_hits'] += 1
            else:
                metadata['ANATOMIA_hits'] += 1
    
    return metadata

def cluster_similar_sequences(sequences, min_length=8, max_length=20):
    """
    Remove redundant overlapping sequences.
    Keep the longest representative from each cluster.
    """
    # Filter by length first
    filtered = {seq: True for seq in sequences if min_length <= len(seq) <= max_length}
    
    # Group overlapping sequences
    clusters = []
    processed = set()
    
    for seq in sorted(filtered.keys(), key=len, reverse=True):
        if seq in processed:
            continue
        
        # Find all sequences that are substrings of this one
        cluster = [seq]
        for other_seq in filtered.keys():
            if other_seq != seq and other_seq not in processed:
                if seq in other_seq or other_seq in seq:
                    cluster.append(other_seq)
                    processed.add(other_seq)
        
        clusters.append(cluster)
        processed.add(seq)
    
    # Return longest from each cluster
    representatives = [max(cluster, key=len) for cluster in clusters]
    
    return representatives

def prioritize_epitopes(output_dir):
    """
    Create prioritized epitope lists with comprehensive metadata.
    """
    
    print("="*60)
    print("EPITOPE PRIORITIZATION & FILTERING")
    print("="*60)
    
    # Load all data
    print("\nLoading data...")
    all_epitopes = pd.read_csv(f"{output_dir}/all_epitopes_with_human_similarity.txt", sep='\t')
    
    # Load source data for metadata
    iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
    
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    exact_df = pd.read_csv(exact_file, sep='\t')
    fuzzy_detailed_df = pd.read_csv(fuzzy_file)
    
    print(f"Total epitopes: {len(all_epitopes)}")
    
    # Filter 1: Human similarity <65%, length 8-20 AA
    print("\nApplying filters...")
    print("  - Human identity <65%")
    print("  - Alignment coverage <70% OR identity <65%")
    print("  - Length 8-20 AA")
    
    filtered = all_epitopes[
        (all_epitopes['Length'] >= 8) & 
        (all_epitopes['Length'] <= 20) &
        (
            (all_epitopes['Max_Human_Identity'] < 65) |
            (all_epitopes['Max_Human_Coverage'] < 70) |
            (all_epitopes['Best_Human_Match'] == 'No_hits')
        )
    ]
    
    print(f"After filtering: {len(filtered)} epitopes")
    
    # Cluster similar/overlapping sequences
    print("\nClustering overlapping sequences...")
    sequences_to_process = filtered['Sequence'].tolist()
    representative_sequences = cluster_similar_sequences(sequences_to_process, 8, 20)
    
    filtered = filtered[filtered['Sequence'].isin(representative_sequences)]
    print(f"After removing overlaps: {len(filtered)} unique epitopes")
    
    # Add comprehensive metadata
    print("\nAdding metadata (this may take a few minutes)...")
    metadata_list = []
    
    for idx, row in filtered.iterrows():
        seq = row['Sequence']
        meta = get_sequence_metadata(seq, iedb_df, exact_df, fuzzy_detailed_df)
        
        metadata_list.append({
            'Epitope_ID': row['Epitope_ID'],
            'Sequence': seq,
            'Length': row['Length'],
            'Source_Methods': row['Source_Methods'],
            'ANATOMIA_Hits': meta['ANATOMIA_hits'],
            'LITERATURE_Hits': meta['LITERATURE_hits'],
            'Total_Hits': meta['ANATOMIA_hits'] + meta['LITERATURE_hits'],
            'Total_Proteins': meta['Total_proteins'],
            'Num_Orthogroups': len(meta['Orthogroups']),
            'Num_Species': len(meta['Species']),
            'Num_IEDB_IDs': len(meta['IEDB_IDs']),
            'Best_Human_Match': row['Best_Human_Match'],
            'Max_Human_Identity': row['Max_Human_Identity'],
            'Max_Human_Coverage': row['Max_Human_Coverage']
        })
    
    enriched_df = pd.DataFrame(metadata_list)
    
    # Calculate priority score
    print("\nCalculating priority scores...")
    
    def calculate_priority(row):
        score = 0
        
        # 1. Total hits (ANATOMIA + LITERATURE) - HIGHEST WEIGHT
        score += row['Total_Hits'] * 100
        
        # 2. Method confidence
        if 'HIGH_CONFIDENCE' in row['Source_Methods']:
            score += 1000
        elif 'Kmer_EXACT' in row['Source_Methods']:
            score += 500
        elif 'IEDB_only' in row['Source_Methods']:
            score += 300
        elif 'Kmer_FUZZY' in row['Source_Methods']:
            score += 100
        
        # 3. Breadth (orthogroups and species)
        score += row['Num_Orthogroups'] * 50
        score += row['Num_Species'] * 30
        
        # 4. Lower human similarity = better
        score += (100 - row['Max_Human_Identity']) * 2
        
        return score
    
    enriched_df['Priority_Score'] = enriched_df.apply(calculate_priority, axis=1)
    enriched_df = enriched_df.sort_values('Priority_Score', ascending=False)
    
    # Create tiered lists
    print("\nCreating tiered epitope lists...")
    
    # Tier 1: HIGH CONFIDENCE + high hits
    tier1 = enriched_df[
        (enriched_df['Source_Methods'].str.contains('HIGH_CONFIDENCE')) &
        (enriched_df['Total_Hits'] >= 5)
    ]
    
    # Tier 2: Kmer_EXACT + moderate hits
    tier2 = enriched_df[
        (enriched_df['Source_Methods'].str.contains('Kmer_EXACT')) &
        (~enriched_df['Source_Methods'].str.contains('HIGH_CONFIDENCE')) &
        (enriched_df['Total_Hits'] >= 3)
    ]
    
    # Tier 3: IEDB or lower hits
    tier3 = enriched_df[
        ~enriched_df['Epitope_ID'].isin(tier1['Epitope_ID']) &
        ~enriched_df['Epitope_ID'].isin(tier2['Epitope_ID'])
    ]
    
    # Save all results
    print("\nSaving results...")
    enriched_df.to_csv(f"{output_dir}/PRIORITIZED_all_epitopes.txt", sep='\t', index=False)
    tier1.to_csv(f"{output_dir}/TIER1_top_priority_epitopes.txt", sep='\t', index=False)
    tier2.to_csv(f"{output_dir}/TIER2_high_priority_epitopes.txt", sep='\t', index=False)
    tier3.to_csv(f"{output_dir}/TIER3_moderate_priority_epitopes.txt", sep='\t', index=False)
    
    # Top 50 overall
    top50 = enriched_df.head(50)
    top50.to_csv(f"{output_dir}/TOP50_epitopes.txt", sep='\t', index=False)
    
    print("\n" + "="*60)
    print("PRIORITIZATION COMPLETE")
    print("="*60)
    print(f"Total filtered epitopes: {len(enriched_df)}")
    print(f"Tier 1 (Top priority): {len(tier1)}")
    print(f"Tier 2 (High priority): {len(tier2)}")
    print(f"Tier 3 (Moderate priority): {len(tier3)}")
    print(f"\nTop 10 epitopes by priority score:")
    print(enriched_df[['Sequence', 'Total_Hits', 'Source_Methods', 'Priority_Score']].head(10).to_string())
    print("="*60)

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
prioritize_epitopes(output_dir)