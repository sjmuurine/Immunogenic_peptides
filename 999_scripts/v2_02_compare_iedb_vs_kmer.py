import pandas as pd
import os
from difflib import SequenceMatcher

def sequence_similarity(seq1, seq2):
    """Calculate similarity between two sequences (0-1)."""
    return SequenceMatcher(None, seq1, seq2).ratio()

def compare_iedb_with_kmers(output_dir):
    """
    Compare IEDB epitopes with k-mer results (both EXACT and FUZZY).
    Checks both high-quality and all IEDB matches.
    """
    
    print("="*80)
    print("COMPARING IEDB EPITOPES WITH K-MER METHODS")
    print("="*80)
    
    # Load k-mer results first (these don't change)
    print("\nLoading k-mer results...")
    
    # K-mer EXACT
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    exact_df = pd.read_csv(exact_file, sep='\t')
    exact_sequences = set(exact_df['Motif'].unique())
    print(f"  K-mer EXACT: {len(exact_sequences)} unique motifs")
    
    # K-mer FUZZY
    fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
    fuzzy_df = pd.read_csv(fuzzy_file)
    fuzzy_sequences = set(fuzzy_df['Consensus'].unique())
    print(f"  K-mer FUZZY: {len(fuzzy_sequences)} unique consensus sequences")
    
    # Function to compare one IEDB file
    def analyze_iedb_file(iedb_file, label):
        """Analyze one IEDB file against k-mer results."""
        
        print("\n" + "="*80)
        print(f"ANALYZING: {label}")
        print("="*80)
        
        iedb_df = pd.read_csv(iedb_file, sep='\t')
        iedb_unique = iedb_df['Epitope_Sequence'].unique()
        
        print(f"\nTotal IEDB records: {len(iedb_df)}")
        print(f"Unique IEDB epitope sequences: {len(iedb_unique)}")
        
        # Check exact matches with k-mer EXACT
        print("\n--- IEDB vs K-mer EXACT (character-by-character match) ---")
        
        exact_matches = []
        for seq in iedb_unique:
            if seq in exact_sequences:
                exact_matches.append(seq)
        
        print(f"IEDB epitopes found in K-mer EXACT: {len(exact_matches)} / {len(iedb_unique)} ({len(exact_matches)/len(iedb_unique)*100:.1f}%)")
        
        # Check matches with k-mer FUZZY (exact match)
        print("\n--- IEDB vs K-mer FUZZY (exact match) ---")
        
        fuzzy_exact_matches = []
        for seq in iedb_unique:
            if seq in fuzzy_sequences:
                fuzzy_exact_matches.append(seq)
        
        print(f"IEDB epitopes found in K-mer FUZZY (exact): {len(fuzzy_exact_matches)} / {len(iedb_unique)} ({len(fuzzy_exact_matches)/len(iedb_unique)*100:.1f}%)")
        
        # Check fuzzy matches (similarity >= 80%)
        print("\n--- IEDB vs K-mer FUZZY (similarity >= 80%) ---")
        print("(This may take a moment...)")
        
        fuzzy_similar_matches = {}
        similarity_threshold = 0.80
        
        for i, iedb_seq in enumerate(iedb_unique):
            if i % 100 == 0:
                print(f"  Checking IEDB sequence {i+1}/{len(iedb_unique)}...", end='\r')
            
            # First check exact match (faster)
            if iedb_seq in fuzzy_sequences:
                fuzzy_similar_matches[iedb_seq] = (iedb_seq, 1.0)
                continue
            
            # Check similar matches (only for same length ±2 for speed)
            best_match = None
            best_similarity = 0
            
            for fuzzy_seq in fuzzy_sequences:
                # Only compare sequences of similar length
                if abs(len(iedb_seq) - len(fuzzy_seq)) <= 2:
                    sim = sequence_similarity(iedb_seq, fuzzy_seq)
                    if sim >= similarity_threshold and sim > best_similarity:
                        best_similarity = sim
                        best_match = fuzzy_seq
            
            if best_match:
                fuzzy_similar_matches[iedb_seq] = (best_match, best_similarity)
        
        print()  # New line after progress
        print(f"IEDB epitopes with fuzzy match (≥80%): {len(fuzzy_similar_matches)} / {len(iedb_unique)} ({len(fuzzy_similar_matches)/len(iedb_unique)*100:.1f}%)")
        
        # Overlap analysis
        print("\n--- METHOD OVERLAP ---")
        
        in_exact_only = set(exact_matches) - set(fuzzy_exact_matches)
        in_fuzzy_only = set(fuzzy_exact_matches) - set(exact_matches)
        in_both = set(exact_matches) & set(fuzzy_exact_matches)
        in_neither = set(iedb_unique) - set(exact_matches) - set(fuzzy_exact_matches)
        
        print(f"IEDB epitopes found in:")
        print(f"  K-mer EXACT only:       {len(in_exact_only):5} ({len(in_exact_only)/len(iedb_unique)*100:5.1f}%)")
        print(f"  K-mer FUZZY only:       {len(in_fuzzy_only):5} ({len(in_fuzzy_only)/len(iedb_unique)*100:5.1f}%)")
        print(f"  BOTH EXACT and FUZZY:   {len(in_both):5} ({len(in_both)/len(iedb_unique)*100:5.1f}%)")
        print(f"  NEITHER (IEDB only):    {len(in_neither):5} ({len(in_neither)/len(iedb_unique)*100:5.1f}%)")
        
        # Create summary dataframe
        summary = []
        for seq in iedb_unique:
            in_exact = seq in exact_matches
            in_fuzzy_exact = seq in fuzzy_exact_matches
            in_fuzzy_similar = seq in fuzzy_similar_matches
            
            fuzzy_match = fuzzy_similar_matches.get(seq, (None, 0))
            
            summary.append({
                'IEDB_Sequence': seq,
                'In_Kmer_EXACT': in_exact,
                'In_Kmer_FUZZY_exact': in_fuzzy_exact,
                'In_Kmer_FUZZY_similar': in_fuzzy_similar,
                'Fuzzy_Best_Match': fuzzy_match[0] if fuzzy_match[0] else 'None',
                'Fuzzy_Similarity': fuzzy_match[1] if fuzzy_match[0] else 0,
                'Method_Category': 'EXACT+FUZZY' if in_exact and in_fuzzy_exact else
                                   'EXACT_only' if in_exact else
                                   'FUZZY_only' if in_fuzzy_exact else
                                   'IEDB_only'
            })
        
        summary_df = pd.DataFrame(summary)
        
        # Save results
        output_file = f"{output_dir}/{label.replace(' ', '_').replace(':', '')}_vs_kmer_comparison.txt"
        summary_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Saved: {output_file}")
        
        return summary_df
    
    # Analyze BOTH IEDB files
    all_iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
    hq_iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders_high_quality.txt'
    
    summary_all = analyze_iedb_file(all_iedb_file, "IEDB_ALL_matches")
    summary_hq = analyze_iedb_file(hq_iedb_file, "IEDB_HIGH_QUALITY_matches")
    
    # Final comparison
    print("\n" + "="*80)
    print("COMPARISON: ALL vs HIGH-QUALITY IEDB")
    print("="*80)
    
    print(f"\nTotal unique IEDB sequences:")
    print(f"  ALL matches:         {len(summary_all)}")
    print(f"  HIGH-QUALITY only:   {len(summary_hq)}")
    print(f"  Difference:          {len(summary_all) - len(summary_hq)}")
    
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    print("\nFor publication, report:")
    print("  1. Use HIGH-QUALITY IEDB matches for downstream analysis")
    print("  2. Report ALL matches for completeness in methods/supplementary")
    print("\n" + "="*80)

# Run analysis
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
os.makedirs(output_dir, exist_ok=True)

compare_iedb_with_kmers(output_dir)