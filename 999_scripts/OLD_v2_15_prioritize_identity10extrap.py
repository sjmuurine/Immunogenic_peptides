import pandas as pd
import sys
import os
from datetime import datetime

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

def is_low_complexity(sequence, threshold=0.4):
    """
    Check if sequence is low complexity (repetitive).
    Returns True if any single amino acid makes up >40% of the sequence.
    """
    if len(sequence) == 0:
        return False
    
    from collections import Counter
    aa_counts = Counter(sequence)
    max_freq = max(aa_counts.values()) / len(sequence)
    
    return max_freq > threshold

def calculate_priority_score(row):
    """
    Calculate priority score for epitope cluster.
    
    Score = (Species × 55) + (IEDB × 20) - (Identity × 30.0) - (Coverage × 0.2)
    
    Higher is better.
    VERY STRONG penalty for human identity - safety prioritized!
    """
    species_score = row['Num_Species'] * 55  # Reduced from 100
    iedb_bonus = 20 if row['Has_IEDB'] else 0
    identity_penalty = row['Max_Human_Identity'] * 30.0  # VERY STRONG PENALTY
    coverage_penalty = row['Max_Human_Coverage'] * 0.2
    
    # Extra penalty for perfect match to human (100% identity)
    perfect_match_penalty = 750 if row['Max_Human_Identity'] == 100.0 else 0
    
    total = species_score + iedb_bonus - identity_penalty - coverage_penalty - perfect_match_penalty
   
    return total

def prioritize_clusters(input_file, output_dir, log_dir):
    """
    Calculate priority scores and filter low-complexity sequences.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_15_prioritize_identity10")
    sys.stdout = logger
    
    print("="*80)
    print("EPITOPE CLUSTER PRIORITIZATION (IDENTITY × 10.0, SPECIES × 90)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load clusters
    print(f"\nLoading clusters: {input_file}")
    clusters_df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total clusters: {len(clusters_df):,}")
    
    # Filter low-complexity sequences
    print("\n" + "="*80)
    print("FILTERING LOW-COMPLEXITY SEQUENCES")
    print("="*80)
    
    clusters_df['Is_Low_Complexity'] = clusters_df['Representative_Sequence'].apply(is_low_complexity)
    
    low_complexity = clusters_df[clusters_df['Is_Low_Complexity'] == True]
    print(f"\nLow-complexity sequences found: {len(low_complexity)}")
    
    if len(low_complexity) > 0:
        print("\nExamples of low-complexity sequences (first 10):")
        for _, row in low_complexity.head(10).iterrows():
            print(f"  {row['Representative_Sequence']:30} - {row['Num_Species']} species")
    
    # Filter them out
    clusters_filtered = clusters_df[clusters_df['Is_Low_Complexity'] == False].copy()
    
    print(f"\nAfter filtering: {len(clusters_filtered):,} clusters")
    print(f"Removed: {len(low_complexity):,} low-complexity sequences")
    
    # Calculate priority scores
    print("\n" + "="*80)
    print("CALCULATING PRIORITY SCORES")
    print("="*80)
    
    print("\nPriority formula:")
    print("  Score = (Species × 90) + (IEDB × 30) - (Identity × 10.0) - (Coverage × 0.5)")
    print("  Where:")
    print("    - Species: Number of species (0-20)")
    print("    - IEDB: 30 if Has_IEDB, 0 otherwise")
    print("    - Identity: Max_Human_Identity (0-100%) - VERY STRONG PENALTY")
    print("    - Coverage: Max_Human_Coverage (0-100%+)")
    print("\n  This formula prioritizes SAFETY (low human similarity) over species breadth")
    
    clusters_filtered['Priority_Score'] = clusters_filtered.apply(calculate_priority_score, axis=1)
    
    # Sort by priority
    clusters_filtered = clusters_filtered.sort_values('Priority_Score', ascending=False)
    
    print(f"\nPriority score statistics:")
    print(f"  Mean: {clusters_filtered['Priority_Score'].mean():.1f}")
    print(f"  Median: {clusters_filtered['Priority_Score'].median():.1f}")
    print(f"  Min: {clusters_filtered['Priority_Score'].min():.1f}")
    print(f"  Max: {clusters_filtered['Priority_Score'].max():.1f}")
    
    # Show top 50
    print("\n" + "="*80)
    print("TOP 50 EPITOPE CLUSTERS")
    print("="*80)
    
    top50 = clusters_filtered.head(50)
    
    print(f"\n{'Rank':<6} {'Sequence':<30} {'Sp':>3} {'IEDB':>4} {'Ident':>6} {'Cov':>6} {'Score':>7}")
    print("-" * 85)
    
    for rank, (_, row) in enumerate(top50.iterrows(), 1):
        seq = row['Representative_Sequence'][:28]
        species = row['Num_Species']
        iedb = 'Yes' if row['Has_IEDB'] else 'No'
        identity = row['Max_Human_Identity']
        coverage = row['Max_Human_Coverage']
        score = row['Priority_Score']
        
        print(f"{rank:<6} {seq:<30} {species:>3} {iedb:>4} {identity:>6.1f} {coverage:>6.1f} {score:>7.1f}")
    
    # Check where TIRIPVHMIETINKLNRISRQML is
    print("\n" + "="*80)
    print("CHECKING TARGET SEQUENCE")
    print("="*80)
    
    target_seq = "TIRIPVHMIETINKLNRISRQML"
    target_row = clusters_filtered[clusters_filtered['Representative_Sequence'] == target_seq]
    
    if len(target_row) > 0:
        rank_idx = target_row.index[0]
        position = clusters_filtered.index.get_loc(rank_idx) + 1
        row = target_row.iloc[0]
        
        print(f"\nSequence: {target_seq}")
        print(f"  Rank: {position}")
        print(f"  Species: {row['Num_Species']}")
        print(f"  IEDB: {row['Has_IEDB']}")
        print(f"  Identity: {row['Max_Human_Identity']:.1f}%")
        print(f"  Coverage: {row['Max_Human_Coverage']:.1f}%")
        print(f"  Score: {row['Priority_Score']:.1f}")
        
        if position <= 50:
            print(f"  ✓ IN TIER1 (TOP 50)!")
        elif position <= 200:
            print(f"  In TIER2 (rank {position})")
        else:
            print(f"  In lower tier (rank {position})")
    else:
        print(f"\nSequence {target_seq} not found in results")
    
    # Create tiers
    print("\n" + "="*80)
    print("CREATING TIERS")
    print("="*80)
    
    top_50 = clusters_filtered.head(50)
    top_200 = clusters_filtered.head(200)
    top_500 = clusters_filtered.head(500)
    top_1000 = clusters_filtered.head(1000)
    
    clusters_filtered['Tier'] = 'TIER5'
    clusters_filtered.loc[clusters_filtered.index.isin(top_1000.index), 'Tier'] = 'TIER4'
    clusters_filtered.loc[clusters_filtered.index.isin(top_500.index), 'Tier'] = 'TIER3'
    clusters_filtered.loc[clusters_filtered.index.isin(top_200.index), 'Tier'] = 'TIER2'
    clusters_filtered.loc[clusters_filtered.index.isin(top_50.index), 'Tier'] = 'TIER1'
    
    print(f"\nTier distribution:")
    tier_counts = clusters_filtered['Tier'].value_counts()
    for tier in ['TIER1', 'TIER2', 'TIER3', 'TIER4', 'TIER5']:
        count = tier_counts.get(tier, 0)
        pct = (count / len(clusters_filtered)) * 100
        print(f"  {tier:6}: {count:5,} clusters ({pct:5.1f}%)")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    all_output = f"{output_dir}/epitope_clusters_prioritized.txt"
    clusters_filtered.to_csv(all_output, sep='\t', index=False)
    print(f"\n✓ Saved all prioritized: {all_output}")
    
    for tier in ['TIER1', 'TIER2', 'TIER3', 'TIER4', 'TIER5']:
        tier_df = clusters_filtered[clusters_filtered['Tier'] == tier]
        if len(tier_df) > 0:
            tier_file = f"{output_dir}/{tier}_clusters.txt"
            tier_df.to_csv(tier_file, sep='\t', index=False)
            print(f"✓ Saved {tier}: {tier_file} ({len(tier_df)} clusters)")
    
    if len(low_complexity) > 0:
        low_complex_file = f"{output_dir}/low_complexity_filtered.txt"
        low_complexity.to_csv(low_complex_file, sep='\t', index=False)
        print(f"✓ Saved low-complexity: {low_complex_file}")
    
    # Analysis
    print("\n" + "="*80)
    print("TIER1 (TOP 50) ANALYSIS")
    print("="*80)
    
    tier1 = clusters_filtered[clusters_filtered['Tier'] == 'TIER1']
    
    print(f"\nSpecies distribution in TIER1:")
    print(f"  Mean: {tier1['Num_Species'].mean():.1f}")
    print(f"  Median: {tier1['Num_Species'].median():.0f}")
    print(f"  Min: {tier1['Num_Species'].min()}")
    print(f"  Max: {tier1['Num_Species'].max()}")
    
    print(f"\nHuman similarity in TIER1:")
    print(f"  Mean identity: {tier1['Max_Human_Identity'].mean():.1f}%")
    print(f"  Median identity: {tier1['Max_Human_Identity'].median():.1f}%")
    print(f"  Min identity: {tier1['Max_Human_Identity'].min():.1f}%")
    print(f"  Max identity: {tier1['Max_Human_Identity'].max():.1f}%")
    
    print(f"\nMethod distribution in TIER1:")
    method_counts = tier1['Representative_Method'].value_counts()
    for method, count in method_counts.items():
        pct = (count / len(tier1)) * 100
        print(f"  {method:30} {count:3} ({pct:5.1f}%)")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_summary.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    prioritize_clusters(input_file, output_dir, log_dir)