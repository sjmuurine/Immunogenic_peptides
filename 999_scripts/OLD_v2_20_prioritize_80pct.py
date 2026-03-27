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

def calculate_priority_score(row):
    """
    Calculate priority score for epitope cluster.
    
    Score = (Species × 100) + (IEDB × 30) - (Identity × 10.0) - (Coverage × 0.2)
    
    Higher is better.
    Balanced formula favoring species breadth and low human identity.
    Should give positive scores since max identity is now <80%.
    """
    species_score = row['Num_Species'] * 100  # Strong species bonus
    iedb_bonus = 30 if row['Has_IEDB'] else 0  # Moderate IEDB bonus
    identity_penalty = row['Max_Human_Identity'] * 10.0  # Strong identity penalty
    coverage_penalty = row['Max_Human_Coverage'] * 0.2  # Light coverage penalty
    
    total = species_score + iedb_bonus - identity_penalty - coverage_penalty
    
    return total

def prioritize_clusters_80pct(input_file, output_dir, log_dir):
    """
    Prioritize high-complexity clusters from 80% identity filtered data.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_20_prioritize_80pct")
    sys.stdout = logger
    
    print("="*80)
    print("EPITOPE CLUSTER PRIORITIZATION (80% FILTERED + HIGH COMPLEXITY)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Load clusters
    print(f"\nLoading: {input_file}")
    clusters_df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total high-complexity clusters: {len(clusters_df):,}")
    
    # Calculate priority scores
    print("\n" + "="*80)
    print("CALCULATING PRIORITY SCORES")
    print("="*80)
    
    print("\nPriority formula:")
    print("  Score = (Species × 100) + (IEDB × 30) - (Identity × 10.0) - (Coverage × 0.2)")
    print("  Where:")
    print("    - Species: Number of species (0-20)")
    print("    - IEDB: 30 if Has_IEDB, 0 otherwise")
    print("    - Identity: Max_Human_Identity (0-80%) - Strong penalty")
    print("    - Coverage: Max_Human_Coverage (0-100%+) - Light penalty")
    
    clusters_df['Priority_Score'] = clusters_df.apply(calculate_priority_score, axis=1)
    
    # Sort by priority
    clusters_df = clusters_df.sort_values('Priority_Score', ascending=False)
    
    print(f"\nPriority score statistics:")
    print(f"  Mean: {clusters_df['Priority_Score'].mean():.1f}")
    print(f"  Median: {clusters_df['Priority_Score'].median():.1f}")
    print(f"  Min: {clusters_df['Priority_Score'].min():.1f}")
    print(f"  Max: {clusters_df['Priority_Score'].max():.1f}")
    
    # Show top 50
    print("\n" + "="*80)
    print("TOP 50 EPITOPE CLUSTERS")
    print("="*80)
    
    top50 = clusters_df.head(50)
    
    print(f"\n{'Rank':<5} {'Sequence':<35} {'Sp':>3} {'ID%':>6} {'IEDB':>4} {'Score':>7}")
    print("-" * 80)
    
    for rank, (_, row) in enumerate(top50.iterrows(), 1):
        seq = row['Representative_Sequence'][:33]
        species = row['Num_Species']
        identity = row['Max_Human_Identity']
        iedb = 'Yes' if row['Has_IEDB'] else 'No'
        score = row['Priority_Score']
        
        print(f"{rank:<5} {seq:<35} {species:>3} {identity:>6.1f} {iedb:>4} {score:>7.1f}")
    
    # Create tiers based on your criteria
    print("\n" + "="*80)
    print("CREATING TIERS")
    print("="*80)
    
    # Your criteria: favor >6 species, but some 3-5 species could be useful
    tier1_candidates = clusters_df[clusters_df['Num_Species'] >= 7]  # Strong candidates
    tier2_candidates = clusters_df[(clusters_df['Num_Species'] >= 5) & (clusters_df['Num_Species'] < 7)]  # Good candidates
    tier3_candidates = clusters_df[(clusters_df['Num_Species'] >= 3) & (clusters_df['Num_Species'] < 5)]  # Potential candidates
    rest_candidates = clusters_df[clusters_df['Num_Species'] < 3]  # Low priority
    
    print(f"\nSpecies-based candidate distribution:")
    print(f"  ≥7 species (TIER1 pool): {len(tier1_candidates):,}")
    print(f"  5-6 species (TIER2 pool): {len(tier2_candidates):,}")
    print(f"  3-4 species (TIER3 pool): {len(tier3_candidates):,}")
    print(f"  <3 species (TIER4 pool): {len(rest_candidates):,}")
    
    # Assign tiers based on priority score within species groups
    clusters_df['Tier'] = 'TIER4'
    
    # TIER1: Top priority from ≥7 species (take top 30)
    if len(tier1_candidates) > 0:
        tier1_top = tier1_candidates.head(30).index
        clusters_df.loc[tier1_top, 'Tier'] = 'TIER1'
    
    # TIER2: Mix of remaining ≥7 species + best 5-6 species (take 50 total)
    tier1_count = (clusters_df['Tier'] == 'TIER1').sum()
    tier2_needed = 50 - tier1_count
    
    tier2_pool = pd.concat([
        tier1_candidates.iloc[30:],  # Remaining ≥7 species
        tier2_candidates  # All 5-6 species
    ]).sort_values('Priority_Score', ascending=False)
    
    if len(tier2_pool) > 0 and tier2_needed > 0:
        tier2_top = tier2_pool.head(tier2_needed).index
        clusters_df.loc[tier2_top, 'Tier'] = 'TIER2'
    
    # TIER3: Best 3-4 species + remaining candidates (take 100 total)
    tier3_assigned = (clusters_df['Tier'].isin(['TIER1', 'TIER2'])).sum()
    tier3_needed = 150 - tier3_assigned
    
    remaining_pool = clusters_df[~clusters_df['Tier'].isin(['TIER1', 'TIER2'])].sort_values('Priority_Score', ascending=False)
    
    if len(remaining_pool) > 0 and tier3_needed > 0:
        tier3_top = remaining_pool.head(tier3_needed).index
        clusters_df.loc[tier3_top, 'Tier'] = 'TIER3'
    
    # Tier distribution
    print(f"\nFinal tier distribution:")
    tier_counts = clusters_df['Tier'].value_counts()
    for tier in ['TIER1', 'TIER2', 'TIER3', 'TIER4']:
        count = tier_counts.get(tier, 0)
        pct = (count / len(clusters_df)) * 100
        print(f"  {tier:6}: {count:5,} clusters ({pct:5.1f}%)")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    # Save all prioritized clusters
    all_output = f"{output_dir}/epitope_clusters_prioritized_80pct.txt"
    clusters_df.to_csv(all_output, sep='\t', index=False)
    print(f"\n✓ Saved all prioritized: {all_output}")
    
    # Save by tier
    for tier in ['TIER1', 'TIER2', 'TIER3']:
        tier_df = clusters_df[clusters_df['Tier'] == tier]
        if len(tier_df) > 0:
            tier_file = f"{output_dir}/{tier}_clusters_80pct.txt"
            tier_df.to_csv(tier_file, sep='\t', index=False)
            print(f"✓ Saved {tier}: {tier_file} ({len(tier_df)} clusters)")
    
    # Analysis
    print("\n" + "="*80)
    print("TIER1 ANALYSIS (TOP CANDIDATES)")
    print("="*80)
    
    tier1 = clusters_df[clusters_df['Tier'] == 'TIER1']
    
    if len(tier1) > 0:
        print(f"\nTIER1 statistics:")
        print(f"  Clusters: {len(tier1)}")
        print(f"  Species range: {tier1['Num_Species'].min()}-{tier1['Num_Species'].max()}")
        print(f"  Mean species: {tier1['Num_Species'].mean():.1f}")
        print(f"  Identity range: {tier1['Max_Human_Identity'].min():.1f}-{tier1['Max_Human_Identity'].max():.1f}%")
        print(f"  Mean identity: {tier1['Max_Human_Identity'].mean():.1f}%")
        print(f"  IEDB sequences: {tier1['Has_IEDB'].sum()}")
        
        print(f"\nMethod distribution in TIER1:")
        method_counts = tier1['Representative_Method'].value_counts()
        for method, count in method_counts.items():
            pct = (count / len(tier1)) * 100
            print(f"  {method:30} {count:3} ({pct:5.1f}%)")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("PRIORITIZATION COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Review TIER1 candidates")
    print("  2. Extract protein regions for top candidates")
    print("  3. Perform alignments and sequence analysis")
    print("  4. Select final epitopes for validation")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_high_complexity.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    prioritize_clusters_80pct(input_file, output_dir, log_dir)