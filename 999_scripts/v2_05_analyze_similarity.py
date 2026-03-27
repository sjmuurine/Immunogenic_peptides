import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def calculate_simple_overlap_similarity(seq1, seq2):
    """Calculate similarity by trying all possible alignments."""
    if seq1 == seq2:
        return 100.0
    
    max_similarity = 0
    
    for offset in range(-len(seq2) + 1, len(seq1)):
        matches = 0
        overlap_length = 0
        
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(seq2):
                overlap_length += 1
                if seq1[i] == seq2[j]:
                    matches += 1
        
        if overlap_length > 0:
            min_len = min(len(seq1), len(seq2))
            if overlap_length >= min_len * 0.8:
                similarity = (matches / overlap_length) * 100
                max_similarity = max(max_similarity, similarity)
    
    return max_similarity

def analyze_similarity_distribution(input_file, output_dir, sample_size=1000):
    """
    Analyze pairwise similarity distribution to choose clustering threshold.
    """
    
    print("="*80)
    print("EPITOPE SIMILARITY DISTRIBUTION ANALYSIS")
    print("="*80)
    
    # Load epitopes
    print(f"\nLoading: {input_file}")
    epitopes_df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total epitopes: {len(epitopes_df)}")
    
    # Filter by length
    epitopes_df = epitopes_df[
        (epitopes_df['Epitope_Sequence'].str.len() >= 11) & 
        (epitopes_df['Epitope_Sequence'].str.len() <= 24)
    ]
    
    print(f"After length filter (11-24 aa): {len(epitopes_df)}")
    
    # Sample for analysis (all pairwise would be too many)
    if len(epitopes_df) > sample_size:
        print(f"\nSampling {sample_size} epitopes for similarity analysis...")
        sample_df = epitopes_df.sample(n=sample_size, random_state=42)
    else:
        sample_df = epitopes_df
    
    sequences = sample_df['Epitope_Sequence'].tolist()
    
    # Calculate pairwise similarities
    print(f"Calculating pairwise similarities for {len(sequences)} sequences...")
    print("(This may take a few minutes...)")
    
    similarities = []
    high_sim_pairs = []
    
    total_pairs = len(sequences) * (len(sequences) - 1) // 2
    calculated = 0
    
    for i in range(len(sequences)):
        if i % 50 == 0:
            print(f"  Progress: {calculated}/{total_pairs} pairs ({calculated/total_pairs*100:.1f}%)", end='\r')
        
        for j in range(i + 1, len(sequences)):
            sim = calculate_simple_overlap_similarity(sequences[i], sequences[j])
            similarities.append(sim)
            calculated += 1
            
            # Track high similarity pairs
            if sim >= 70:
                high_sim_pairs.append({
                    'Seq1': sequences[i],
                    'Seq2': sequences[j],
                    'Similarity': sim
                })
    
    print()
    
    # Statistics
    similarities = np.array(similarities)
    
    print("\n" + "="*80)
    print("SIMILARITY STATISTICS")
    print("="*80)
    
    print(f"\nTotal pairwise comparisons: {len(similarities)}")
    print(f"Mean similarity: {similarities.mean():.1f}%")
    print(f"Median similarity: {np.median(similarities):.1f}%")
    print(f"Std deviation: {similarities.std():.1f}%")
    
    print("\nSimilarity distribution:")
    thresholds = [70, 75, 80, 85, 90, 95, 100]
    for thresh in thresholds:
        count = (similarities >= thresh).sum()
        pct = (count / len(similarities)) * 100
        print(f"  >= {thresh}%: {count:8} pairs ({pct:5.2f}%)")
    
    # Estimate cluster sizes at different thresholds
    print("\n" + "="*80)
    print("ESTIMATED IMPACT OF DIFFERENT THRESHOLDS")
    print("="*80)
    
    print("\nIf we cluster at different thresholds:")
    for thresh in [70, 75, 80, 85, 90]:
        high_sim = (similarities >= thresh).sum()
        # Rough estimate: if X% of pairs are similar, clusters will be ~sqrt(X%) of original
        est_reduction = (high_sim / len(similarities)) ** 0.5
        est_clusters = int(len(epitopes_df) * (1 - est_reduction * 0.5))
        
        print(f"  {thresh}% threshold:")
        print(f"    High-similarity pairs: {high_sim} ({high_sim/len(similarities)*100:.2f}%)")
        print(f"    Estimated clusters: ~{est_clusters:,}")
    
    # Show examples of high-similarity pairs
    print("\n" + "="*80)
    print("EXAMPLES OF HIGH-SIMILARITY PAIRS")
    print("="*80)
    
    high_sim_df = pd.DataFrame(high_sim_pairs)
    
    if len(high_sim_df) > 0:
        # Group by similarity ranges
        for lower, upper in [(95, 100), (90, 95), (85, 90), (80, 85), (75, 80), (70, 75)]:
            range_pairs = high_sim_df[
                (high_sim_df['Similarity'] >= lower) & 
                (high_sim_df['Similarity'] < upper)
            ]
            
            if len(range_pairs) > 0:
                print(f"\n{lower}-{upper}% similarity ({len(range_pairs)} pairs):")
                for _, row in range_pairs.head(3).iterrows():
                    print(f"  {row['Similarity']:.1f}%:")
                    print(f"    {row['Seq1']}")
                    print(f"    {row['Seq2']}")
    
    # Create visualization
    print("\n" + "="*80)
    print("CREATING VISUALIZATION")
    print("="*80)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax1 = axes[0]
    ax1.hist(similarities, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Pairwise Similarity (%)', fontsize=12)
    ax1.set_ylabel('Number of Pairs', fontsize=12)
    ax1.set_title('Distribution of Pairwise Similarities', fontsize=14, fontweight='bold')
    ax1.axvline(x=70, color='red', linestyle='--', label='70% threshold')
    ax1.axvline(x=80, color='orange', linestyle='--', label='80% threshold')
    ax1.axvline(x=90, color='green', linestyle='--', label='90% threshold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Cumulative
    ax2 = axes[1]
    sorted_sim = np.sort(similarities)
    cumulative = np.arange(1, len(sorted_sim) + 1) / len(sorted_sim) * 100
    ax2.plot(sorted_sim, cumulative, linewidth=2, color='steelblue')
    ax2.set_xlabel('Similarity Threshold (%)', fontsize=12)
    ax2.set_ylabel('Cumulative % of Pairs Above Threshold', fontsize=12)
    ax2.set_title('Cumulative Distribution', fontsize=14, fontweight='bold')
    ax2.axvline(x=70, color='red', linestyle='--', alpha=0.5)
    ax2.axvline(x=80, color='orange', linestyle='--', alpha=0.5)
    ax2.axvline(x=90, color='green', linestyle='--', alpha=0.5)
    ax2.grid(alpha=0.3)
    ax2.invert_yaxis()
    
    plt.tight_layout()
    
    plot_file = f"{output_dir}/similarity_distribution_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {plot_file}")
    plt.close()
    
    # Save high-similarity pairs
    if len(high_sim_df) > 0:
        pairs_file = f"{output_dir}/high_similarity_pairs.txt"
        high_sim_df.sort_values('Similarity', ascending=False).to_csv(pairs_file, sep='\t', index=False)
        print(f"✓ Saved: {pairs_file}")
    
    print("\n" + "="*80)
    print("RECOMMENDATION")
    print("="*80)
    
    pct_90 = (similarities >= 90).sum() / len(similarities) * 100
    pct_80 = (similarities >= 80).sum() / len(similarities) * 100
    pct_70 = (similarities >= 70).sum() / len(similarities) * 100
    
    print(f"\nBased on the analysis:")
    print(f"  - 90% threshold: {pct_90:.2f}% of pairs cluster (strictest)")
    print(f"  - 80% threshold: {pct_80:.2f}% of pairs cluster (moderate)")
    print(f"  - 70% threshold: {pct_70:.2f}% of pairs cluster (permissive)")
    
    if pct_90 < 1:
        print("\n→ 90% threshold may be too strict (very few clusters formed)")
        print("→ Consider 75-80% threshold")
    elif pct_70 > 10:
        print("\n→ 70% threshold may be too permissive (merges distant sequences)")
        print("→ Consider 80-85% threshold")
    else:
        print("\n→ 75-85% threshold range seems reasonable")
    
    print("="*80)

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    
    # Analyze with 1000 random epitopes
    analyze_similarity_distribution(input_file, output_dir, sample_size=1000)