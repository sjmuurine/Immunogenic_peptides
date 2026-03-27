import pandas as pd
import re
from collections import Counter

def calculate_sequence_complexity(seq):
    """
    Calculate sequence complexity (0-1, higher = more complex).
    Low complexity = repetitive sequences.
    """
    if len(seq) == 0:
        return 0
    
    # Count unique amino acids
    unique_aa = len(set(seq))
    
    # Penalize if too few unique AAs
    if unique_aa < 5:
        base_complexity = unique_aa / 10
    else:
        base_complexity = 0.5 + (unique_aa / 40)
    
    # Check for obvious repeats
    has_dinucleotide_repeat = bool(re.search(r'(.{2})\1{3,}', seq))  # AA repeated 4+ times
    has_single_repeat = bool(re.search(r'(.)\1{5,}', seq))           # Single AA 6+ times
    
    # Penalize repeats heavily
    if has_dinucleotide_repeat:
        base_complexity *= 0.3
    if has_single_repeat:
        base_complexity *= 0.3
    
    return base_complexity

def is_low_complexity(seq, threshold=0.25):
    """
    Check if sequence is low complexity (repetitive).
    Lower threshold = more permissive (keep more sequences)
    """
    complexity = calculate_sequence_complexity(seq)
    
    # Also check for specific problematic patterns
    # Too many of same amino acid
    aa_counts = Counter(seq)
    most_common_aa, most_common_count = aa_counts.most_common(1)[0]
    if most_common_count / len(seq) > 0.6:  # >60% single amino acid
        return True
    
    return complexity < threshold

def recalculate_priority_score(row):
    """
    Better priority scoring that emphasizes breadth over raw counts.
    """
    score = 0
    
    # 1. Species and Orthogroup diversity - HIGHEST WEIGHT
    score += row['Unique_Species_AllVariants'] * 200
    score += row['Unique_OGs_AllVariants'] * 150
    
    # 2. Method confidence
    if row['Has_HIGH_CONFIDENCE']:
        score += 2000
    
    if row['Has_IEDB']:
        score += 800
    
    if row['Has_Kmer_EXACT']:
        score += 500
    else:
        score += 200
    
    # 3. Total hits - moderate weight
    score += row['Total_Hits_AllVariants'] * 30
    
    # 4. Both ANATOMIA and LITERATURE
    has_both = (row['ANATOMIA_Hits_AllVariants'] > 0 and 
                row['LITERATURE_Hits_AllVariants'] > 0)
    if has_both:
        score += 1000
    
    # 5. Total proteins
    score += row['Total_Proteins_AllVariants'] * 20
    
    # 6. Lower human similarity
    score += (100 - row['Human_Identity']) * 3
    
    # 7. PENALTY for poor breadth despite high hits
    if row['Total_Hits_AllVariants'] > 20 and row['Unique_Species_AllVariants'] <= 1:
        score -= 2000
    
    # 8. Complexity penalty
    complexity = calculate_sequence_complexity(row['Representative_Sequence'])
    if complexity < 0.25:
        score *= 0.5
    
    return score

def rerank_epitope_clusters(cluster_dir, output_dir):
    """
    Re-rank clusters with improved scoring and complexity filtering.
    """
    
    print("="*70)
    print("RE-RANKING EPITOPE CLUSTERS")
    print("="*70)
    
    datasets = ['TOP50', 'TIER1', 'TIER2', 'COMBINED_TOP_TIERS']
    
    for dataset in datasets:
        cluster_file = f"{cluster_dir}/{dataset}_cluster_summary_v2.txt"
        
        try:
            df = pd.read_csv(cluster_file, sep='\t')
        except:
            print(f"\nSkipping {dataset} (file not found)")
            continue
        
        print(f"\n{'='*70}")
        print(f"Processing: {dataset}")
        print(f"{'='*70}")
        print(f"Total clusters: {len(df)}")
        
        # Calculate complexity
        df['Sequence_Complexity'] = df['Representative_Sequence'].apply(
            calculate_sequence_complexity
        )
        
        df['Is_Low_Complexity'] = df['Representative_Sequence'].apply(
            lambda x: is_low_complexity(x, threshold=0.25)
        )
        
        # Filter
        high_complexity = df[~df['Is_Low_Complexity']].copy()
        low_complexity = df[df['Is_Low_Complexity']].copy()
        
        print(f"High complexity clusters: {len(high_complexity)}")
        print(f"Low complexity (filtered): {len(low_complexity)}")
        
        if len(high_complexity) == 0:
            print("WARNING: All sequences filtered! Using more permissive threshold...")
            # Use top 80% by complexity
            threshold_80 = df['Sequence_Complexity'].quantile(0.2)
            high_complexity = df[df['Sequence_Complexity'] >= threshold_80].copy()
            low_complexity = df[df['Sequence_Complexity'] < threshold_80].copy()
            print(f"Adjusted - High complexity: {len(high_complexity)}, Low: {len(low_complexity)}")
        
        # Recalculate priority scores (fixed - return single value)
        print("Calculating new priority scores...")
        new_scores = []
        for idx, row in high_complexity.iterrows():
            score = recalculate_priority_score(row)
            new_scores.append(score)
        
        high_complexity['New_Priority_Score'] = new_scores
        
        # Sort
        high_complexity = high_complexity.sort_values('New_Priority_Score', ascending=False)
        
        # Save
        output_file = f"{output_dir}/{dataset}_RERANKED.txt"
        high_complexity.to_csv(output_file, sep='\t', index=False)
        
        if len(low_complexity) > 0:
            filtered_file = f"{output_dir}/{dataset}_filtered_low_complexity.txt"
            low_complexity.to_csv(filtered_file, sep='\t', index=False)
        
        print(f"\n✓ Saved re-ranked: {output_file}")
        if len(low_complexity) > 0:
            print(f"✓ Saved filtered: {filtered_file}")
        
        # Show top 10
        print(f"\nTop 10 after re-ranking:")
        top10 = high_complexity.head(10)[
            ['Cluster_ID', 'Representative_Sequence', 'Cluster_Size',
             'Unique_Species_AllVariants', 'Unique_OGs_AllVariants', 
             'Total_Hits_AllVariants', 'Methods_Used', 'New_Priority_Score']
        ]
        print(top10.to_string(index=False))
        
        # Show what got filtered
        if len(low_complexity) > 0:
            print(f"\nExamples of filtered low-complexity sequences:")
            for seq in low_complexity.head(5)['Representative_Sequence']:
                complexity = calculate_sequence_complexity(seq)
                print(f"  {seq} (complexity: {complexity:.2f})")
    
    print("\n" + "="*70)
    print("RE-RANKING COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {output_dir}")
    print("="*70)

# Run
cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
rerank_epitope_clusters(cluster_dir, output_dir)