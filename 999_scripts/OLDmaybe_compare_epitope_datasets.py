import pandas as pd
import os

def compare_epitope_datasets(iedb_file, exact_file, fuzzy_file, output_dir):
    """
    Compare IEDB-matched epitopes with k-mer predicted epitopes.
    Finds overlaps and unique sequences from each method.
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("EPITOPE DATASET COMPARISON")
    print("="*60)
    
    # Load IEDB results
    print("\nLoading IEDB BLAST results...")
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    iedb_sequences = set(iedb_df['Epitope_Sequence'].unique())
    print(f"  Unique IEDB epitope sequences: {len(iedb_sequences)}")
    
    # Load exact match k-mer results
    print("\nLoading exact match k-mer results...")
    exact_df = pd.read_csv(exact_file, sep='\t')
    exact_sequences = set(exact_df['Motif'].unique())
    print(f"  Unique exact match motifs: {len(exact_sequences)}")
    
    # Load fuzzy match k-mer results (use SUMMARY file, not detailed file)
    print("\nLoading fuzzy match k-mer results...")
    fuzzy_summary_file = fuzzy_file.replace('epitope_detailed.csv', 'epitope_summary.csv')
    fuzzy_summary_df = pd.read_csv(fuzzy_summary_file)
    fuzzy_detailed_df = pd.read_csv(fuzzy_file)  # Keep this for detailed info if needed
    fuzzy_sequences = set(fuzzy_summary_df['Consensus_Sequence'].unique())
    print(f"  Unique fuzzy consensus sequences: {len(fuzzy_sequences)}")
    
    # Combine k-mer results (remove duplicates)
    kmer_sequences = exact_sequences.union(fuzzy_sequences)
    print(f"\nTotal unique k-mer sequences (exact + fuzzy): {len(kmer_sequences)}")
    
    # Find overlaps
    iedb_and_exact = iedb_sequences.intersection(exact_sequences)
    iedb_and_fuzzy = iedb_sequences.intersection(fuzzy_sequences)
    iedb_and_any_kmer = iedb_sequences.intersection(kmer_sequences)
    
    print("\n" + "="*60)
    print("OVERLAP ANALYSIS")
    print("="*60)
    print(f"IEDB sequences also in EXACT k-mer: {len(iedb_and_exact)}")
    print(f"IEDB sequences also in FUZZY k-mer: {len(iedb_and_fuzzy)}")
    print(f"IEDB sequences in ANY k-mer method: {len(iedb_and_any_kmer)}")
    print(f"IEDB sequences ONLY (not in k-mer): {len(iedb_sequences - kmer_sequences)}")
    print(f"K-mer sequences ONLY (not in IEDB): {len(kmer_sequences - iedb_sequences)}")
    
    # Create detailed comparison tables
    
    # 1. Sequences found by both methods (HIGH CONFIDENCE)
    overlap_data = []
    for seq in iedb_and_any_kmer:
        in_exact = seq in exact_sequences
        in_fuzzy = seq in fuzzy_sequences
        
        # Get IEDB info
        iedb_info = iedb_df[iedb_df['Epitope_Sequence'] == seq]
        iedb_ids = '; '.join(iedb_info['Epitope_ID'].unique()[:3])  # First 3
        anatomia_hits = len(iedb_info[iedb_info['Data_Source'] == 'ANATOMIA'])
        literature_hits = len(iedb_info[iedb_info['Data_Source'] == 'LITERATURE'])
        
        # Get k-mer info
        if in_exact:
            exact_info = exact_df[exact_df['Motif'] == seq].iloc[0]
            kmer_count = exact_info['UniqueSeqCount']
            fold_enrich = exact_info['FoldEnrichment']
        elif in_fuzzy:
            fuzzy_info = fuzzy_summary_df[fuzzy_summary_df['Consensus_Sequence'] == seq].iloc[0]
            kmer_count = fuzzy_info['Occurrences']
            fold_enrich = "N/A (fuzzy)"
        
        overlap_data.append({
            'Sequence': seq,
            'Length': len(seq),
            'Found_In_Methods': 'IEDB+Exact' if in_exact else 'IEDB+Fuzzy',
            'IEDB_IDs': iedb_ids,
            'IEDB_ANATOMIA_Hits': anatomia_hits,
            'IEDB_LITERATURE_Hits': literature_hits,
            'Kmer_Occurrences': kmer_count,
            'Kmer_Fold_Enrichment': fold_enrich
        })
    
    overlap_df = pd.DataFrame(overlap_data)
    overlap_df = overlap_df.sort_values('IEDB_ANATOMIA_Hits', ascending=False)
    overlap_df.to_csv(f"{output_dir}/HIGH_CONFIDENCE_overlapping_epitopes.txt", sep='\t', index=False)
    print(f"\n✓ High confidence overlaps saved: HIGH_CONFIDENCE_overlapping_epitopes.txt")
    
    # 2. IEDB-only sequences
    iedb_only = []
    for seq in (iedb_sequences - kmer_sequences):
        iedb_info = iedb_df[iedb_df['Epitope_Sequence'] == seq]
        iedb_ids = '; '.join(iedb_info['Epitope_ID'].unique()[:3])
        anatomia_hits = len(iedb_info[iedb_info['Data_Source'] == 'ANATOMIA'])
        literature_hits = len(iedb_info[iedb_info['Data_Source'] == 'LITERATURE'])
        
        iedb_only.append({
            'Sequence': seq,
            'Length': len(seq),
            'IEDB_IDs': iedb_ids,
            'ANATOMIA_Hits': anatomia_hits,
            'LITERATURE_Hits': literature_hits
        })
    
    iedb_only_df = pd.DataFrame(iedb_only)
    iedb_only_df = iedb_only_df.sort_values('ANATOMIA_Hits', ascending=False)
    iedb_only_df.to_csv(f"{output_dir}/IEDB_only_epitopes.txt", sep='\t', index=False)
    print(f"✓ IEDB-only sequences saved: IEDB_only_epitopes.txt")
    
    # 3. K-mer only sequences (from exact matches)
    kmer_only_exact = []
    for seq in (exact_sequences - iedb_sequences):
        exact_info = exact_df[exact_df['Motif'] == seq].iloc[0]
        kmer_only_exact.append({
            'Sequence': seq,
            'Length': len(seq),
            'UniqueSeqCount': exact_info['UniqueSeqCount'],
            'FoldEnrichment': exact_info['FoldEnrichment'],
            'Source_Proteins': str(exact_info['Sequence_IDs'])[:200]  # Truncate if too long
        })
    
    kmer_only_exact_df = pd.DataFrame(kmer_only_exact)
    kmer_only_exact_df = kmer_only_exact_df.sort_values('UniqueSeqCount', ascending=False)
    kmer_only_exact_df.to_csv(f"{output_dir}/Kmer_EXACT_only_epitopes.txt", sep='\t', index=False)
    print(f"✓ K-mer exact-only sequences saved: Kmer_EXACT_only_epitopes.txt")
    
    # 4. K-mer fuzzy only
    kmer_only_fuzzy = []
    for seq in (fuzzy_sequences - iedb_sequences - exact_sequences):
        fuzzy_info = fuzzy_summary_df[fuzzy_summary_df['Consensus_Sequence'] == seq].iloc[0]
        kmer_only_fuzzy.append({
            'Sequence': seq,
            'Length': fuzzy_info['Length'],
            'Occurrences': fuzzy_info['Occurrences'],
            'Unique_Sequences': fuzzy_info['Unique_Sequences'],
            'Avg_Similarity': fuzzy_info['Avg_Similarity']
        })
    
    kmer_only_fuzzy_df = pd.DataFrame(kmer_only_fuzzy)
    kmer_only_fuzzy_df = kmer_only_fuzzy_df.sort_values('Occurrences', ascending=False)
    kmer_only_fuzzy_df.to_csv(f"{output_dir}/Kmer_FUZZY_only_epitopes.txt", sep='\t', index=False)
    print(f"✓ K-mer fuzzy-only sequences saved: Kmer_FUZZY_only_epitopes.txt")
    
    # Summary statistics
    summary = {
        'Category': [
            'Total IEDB epitopes',
            'Total K-mer epitopes (exact)',
            'Total K-mer epitopes (fuzzy)',
            'Total K-mer epitopes (combined)',
            'High confidence (IEDB + K-mer)',
            'IEDB only',
            'K-mer only (exact)',
            'K-mer only (fuzzy)',
            'All unique epitopes'
        ],
        'Count': [
            len(iedb_sequences),
            len(exact_sequences),
            len(fuzzy_sequences),
            len(kmer_sequences),
            len(iedb_and_any_kmer),
            len(iedb_sequences - kmer_sequences),
            len(exact_sequences - iedb_sequences),
            len(fuzzy_sequences - iedb_sequences - exact_sequences),
            len(iedb_sequences.union(kmer_sequences))
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f"{output_dir}/comparison_summary.txt", sep='\t', index=False)
    print(f"✓ Summary saved: comparison_summary.txt")
    
    print("\n" + "="*60)
    print("COMPARISON COMPLETE")
    print("="*60)
    print(f"Results saved to: {output_dir}")
    
    return overlap_df, iedb_only_df, kmer_only_exact_df, kmer_only_fuzzy_df

# ============================================
# Run the comparison
# ============================================

iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'

compare_epitope_datasets(iedb_file, exact_file, fuzzy_file, output_dir)