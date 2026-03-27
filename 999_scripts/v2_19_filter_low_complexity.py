import pandas as pd
import sys
import os
from datetime import datetime
from collections import Counter

class Logger:
    """Dual output to console and log file."""
    def __init__(self, log_dir, script_name):
        os.makedirs(log_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = f"{log_dir}/{timestamp}_{script_name}.log"
        self.terminal = sys.stdout
        self.log = open(self.log_file, 'w')
        print(f"Logging to: {self.log_file}")
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.close()
        sys.stdout = self.terminal

def is_low_complexity(sequence, max_single_aa=0.4, max_dinucleotide=0.6):
    """
    Check if sequence is low complexity.
    
    Returns True if:
    - Any single amino acid makes up >40% of sequence
    - Any dinucleotide repeat makes up >60% of sequence
    """
    if len(sequence) <= 3:
        return True
    
    # Single amino acid frequency
    aa_counts = Counter(sequence)
    max_aa_freq = max(aa_counts.values()) / len(sequence)
    
    if max_aa_freq > max_single_aa:
        return True
    
    # Dinucleotide repeats
    if len(sequence) >= 4:
        dinucs = [sequence[i:i+2] for i in range(len(sequence)-1)]
        dinuc_counts = Counter(dinucs)
        max_dinuc_freq = max(dinuc_counts.values()) / len(dinucs)
        
        if max_dinuc_freq > max_dinucleotide:
            return True
    
    return False

def filter_low_complexity(input_file, output_dir, log_dir):
    """
    Filter out low-complexity sequences from cluster results.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_19_filter_low_complexity")
    sys.stdout = logger
    
    print("="*80)
    print("LOW-COMPLEXITY SEQUENCE FILTER")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nFiltering Criteria:")
    print("  - Single amino acid >40% of sequence")
    print("  - Dinucleotide repeat >60% of sequence")
    print("  - Sequences ≤3 amino acids")
    
    # Load clusters
    print(f"\nLoading: {input_file}")
    clusters_df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total clusters: {len(clusters_df):,}")
    
    # Check complexity
    print("\n" + "="*80)
    print("ANALYZING SEQUENCE COMPLEXITY")
    print("="*80)
    
    clusters_df['Is_Low_Complexity'] = clusters_df['Representative_Sequence'].apply(is_low_complexity)
    
    low_complexity = clusters_df[clusters_df['Is_Low_Complexity'] == True]
    high_complexity = clusters_df[clusters_df['Is_Low_Complexity'] == False]
    
    print(f"\nComplexity analysis:")
    print(f"  Low complexity: {len(low_complexity):,} ({len(low_complexity)/len(clusters_df)*100:.1f}%)")
    print(f"  High complexity: {len(high_complexity):,} ({len(high_complexity)/len(clusters_df)*100:.1f}%)")
    
    # Show examples of low complexity
    if len(low_complexity) > 0:
        print(f"\nExamples of LOW complexity sequences (filtered out):")
        for _, row in low_complexity.head(10).iterrows():
            seq = row['Representative_Sequence']
            species = row['Num_Species']
            # Calculate dominant AA
            aa_counts = Counter(seq)
            dominant_aa, count = aa_counts.most_common(1)[0]
            freq = count / len(seq) * 100
            print(f"  {seq:25} ({species}sp) - {dominant_aa}: {freq:.0f}%")
    
    # Statistics on high complexity sequences
    print(f"\n" + "="*80)
    print("HIGH COMPLEXITY SEQUENCES (KEPT)")
    print("="*80)
    
    print(f"\nSpecies distribution:")
    print(f"  Mean: {high_complexity['Num_Species'].mean():.1f}")
    print(f"  Median: {high_complexity['Num_Species'].median():.0f}")
    print(f"  Max: {high_complexity['Num_Species'].max()}")
    
    print(f"\nHuman identity distribution:")
    print(f"  Mean: {high_complexity['Max_Human_Identity'].mean():.1f}%")
    print(f"  Median: {high_complexity['Max_Human_Identity'].median():.1f}%")
    print(f"  Min: {high_complexity['Max_Human_Identity'].min():.1f}%")
    print(f"  Max: {high_complexity['Max_Human_Identity'].max():.1f}%")
    
    # Show identity distribution
    for threshold in [0, 20, 40, 60, 70, 75]:
        count = (high_complexity['Max_Human_Identity'] >= threshold).sum()
        pct = count / len(high_complexity) * 100
        print(f"  ≥{threshold}%: {count:,} ({pct:.1f}%)")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    # Save high complexity (kept)
    high_complexity_file = f"{output_dir}/epitope_clusters_high_complexity.txt"
    high_complexity.to_csv(high_complexity_file, sep='\t', index=False)
    print(f"\n✓ Saved HIGH complexity (kept): {high_complexity_file}")
    print(f"  Clusters: {len(high_complexity):,}")
    
    # Save low complexity (filtered out)
    low_complexity_file = f"{output_dir}/epitope_clusters_low_complexity_filtered.txt"
    low_complexity.to_csv(low_complexity_file, sep='\t', index=False)
    print(f"✓ Saved LOW complexity (filtered): {low_complexity_file}")
    print(f"  Clusters: {len(low_complexity):,}")
    
    print(f"\nTop 20 high-complexity clusters:")
    top20 = high_complexity.head(20)
    print(f"\n{'Rank':<4} {'Sequence':<35} {'Sp':<3} {'ID%':<5} {'Method'}")
    print("-" * 70)
    for i, (_, row) in enumerate(top20.iterrows(), 1):
        seq = row['Representative_Sequence'][:33]
        species = row['Num_Species']
        identity = row['Max_Human_Identity']
        method = 'IEDB' if row['Has_IEDB'] else 'Exact'
        print(f"{i:<4} {seq:<35} {species:<3} {identity:<5.1f} {method}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("LOW-COMPLEXITY FILTERING COMPLETE")
    print("="*80)
    
    logger.close()
    
    return high_complexity_file

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_80pct_summary.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    filter_low_complexity(input_file, output_dir, log_dir)