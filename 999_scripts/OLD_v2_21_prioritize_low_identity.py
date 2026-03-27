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
    Calculate priority score emphasizing low human identity.
    
    Score = (Species × 50) + (IEDB × 30) - (Identity × 20.0) + (Low_Identity_Bonus)
    
    Higher is better.
    Strong emphasis on bacteria-specific epitopes with low human similarity.
    """
    species_score = row['Num_Species'] * 50  # Reduced from 100
    iedb_bonus = 30 if row['Has_IEDB'] else 0
    identity_penalty = row['Max_Human_Identity'] * 20.0  # Strong penalty
    
    # Extra bonus for very low identity (bacteria-specific)
    identity = row['Max_Human_Identity']
    low_identity_bonus = 0
    if identity < 30:
        low_identity_bonus = 200  # Big bonus for very bacteria-specific
    elif identity < 50:
        low_identity_bonus = 100  # Good bonus for low identity
    elif identity < 60:
        low_identity_bonus = 50   # Small bonus for moderate identity
    
    total = species_score + iedb_bonus - identity_penalty + low_identity_bonus
    
    return total

def prioritize_clusters_low_identity(input_file, output_dir, log_dir):
    """
    Prioritize clusters with emphasis on low human identity.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_21_prioritize_low_identity")
    sys.stdout = logger
    
    print("="*80)
    print("EPITOPE PRIORITIZATION (LOW HUMAN IDENTITY FOCUS)")
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
    
    print("\nPriority formula (LOW IDENTITY FOCUS):")
    print("  Score = (Species × 50) + (IEDB × 30) - (Identity × 20.0) + Low_Identity_Bonus")
    print("  Where Low_Identity_Bonus:")
    print("    - <30% identity: +200 (very bacteria-specific)")
    print("    - 30-50% identity: +100 (bacteria-specific)")
    print("    - 50-60% identity: +50 (moderately specific)")
    print("    - >60% identity: +0 (no bonus)")
    
    clusters_df['Priority_Score'] = clusters_df.apply(calculate_priority_score, axis=1)
    
    # Sort by priority
    clusters_df = clusters_df.sort_values('Priority_Score', ascending=False)
    
    print(f"\nPriority score statistics:")
    print(f"  Mean: {clusters_df['Priority_Score'].mean():.1f}")
    print(f"  Median: {clusters_df['Priority_Score'].median():.1f}")
    print(f"  Min: {clusters_df['Priority_Score'].min():.1f}")
    print(f"  Max: {clusters_df['Priority_Score'].max():.1f}")
    
    # Show identity distribution impact
    print(f"\nIdentity distribution and bonuses:")
    for threshold, bonus in [(30, 200), (50, 100), (60, 50)]:
        count = (clusters_df['Max_Human_Identity'] < threshold).sum()
        pct = count / len(clusters_df) * 100
        print(f"  <{threshold}% identity (bonus +{bonus}): {count:,} clusters ({pct:.1f}%)")
    
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
    
    # Analyze how the new scoring affects different categories
    print("\n" + "="*80)
    print("SCORING IMPACT ANALYSIS")
    print("="*80)
    
    print(f"\nTop 50 characteristics:")
    print(f"  Species range: {top50['Num_Species'].min()}-{top50['Num_Species'].max()}")
    print(f"  Mean species: {top50['Num_Species'].mean():.1f}")
    print(f"  Identity range: {top50['Max_Human_Identity'].min():.1f}-{top50['Max_Human_Identity'].max():.1f}%")
    print(f"  Mean identity: {top50['Max_Human_Identity'].mean():.1f}%")
    print(f"  IEDB sequences: {top50['Has_IEDB'].sum()}")
    
    print(f"\nSpecies distribution in TOP 50:")
    species_dist = top50['Num_Species'].value_counts().sort_index(ascending=False)
    for species_count, freq in species_dist.items():
        print(f"  {species_count} species: {freq} epitopes")
    
    print(f"\nIdentity distribution in TOP 50:")
    for threshold in [30, 40, 50, 60, 70]:
        count = (top50['Max_Human_Identity'] < threshold).sum()
        pct = count / len(top50) * 100
        print(f"  <{threshold}%: {count} ({pct:.0f}%)")
    
    # Create simplified tiering (2 tiers as requested)
    print("\n" + "="*80)
    print("CREATING SIMPLIFIED TIERS (2 TIERS)")
    print("="*80)
    
    # TIER1: Top 50 for intensive analysis
    tier1_indices = clusters_df.head(50).index
    
    # TIER2: Next 100 as backup candidates
    tier2_indices = clusters_df.iloc[50:150].index
    
    clusters_df['Tier'] = 'TIER3'  # Default
    clusters_df.loc[tier1_indices, 'Tier'] = 'TIER1'
    clusters_df.loc[tier2_indices, 'Tier'] = 'TIER2'
    
    print(f"\nTier distribution:")
    tier_counts = clusters_df['Tier'].value_counts()
    for tier in ['TIER1', 'TIER2', 'TIER3']:
        count = tier_counts.get(tier, 0)
        pct = (count / len(clusters_df)) * 100
        print(f"  {tier:6}: {count:5,} clusters ({pct:5.1f}%)")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    # Save all prioritized clusters
    all_output = f"{output_dir}/epitope_clusters_prioritized_low_identity.txt"
    clusters_df.to_csv(all_output, sep='\t', index=False)
    print(f"\n✓ Saved all prioritized: {all_output}")
    
    # Save by tier
    for tier in ['TIER1', 'TIER2']:
        tier_df = clusters_df[clusters_df['Tier'] == tier]
        if len(tier_df) > 0:
            tier_file = f"{output_dir}/{tier}_clusters_low_identity.txt"
            tier_df.to_csv(tier_file, sep='\t', index=False)
            print(f"✓ Saved {tier}: {tier_file} ({len(tier_df)} clusters)")
    
    # Show some specific examples of interest
    print("\n" + "="*80)
    print("EXAMPLES OF INTEREST")
    print("="*80)
    
    # Look for low-species, low-identity candidates
    interesting = clusters_df[
        (clusters_df['Num_Species'] >= 3) & 
        (clusters_df['Num_Species'] <= 6) & 
        (clusters_df['Max_Human_Identity'] < 50)
    ].head(10)
    
    print(f"\nLow-species, low-identity candidates (bacteria-specific):")
    if len(interesting) > 0:
        print(f"{'Sequence':<35} {'Sp':<3} {'ID%':<6} {'Species'}")
        print("-" * 80)
        for _, row in interesting.iterrows():
            seq = row['Representative_Sequence'][:33]
            species = row['Num_Species']
            identity = row['Max_Human_Identity']
            species_list = row['Species_List'][:40]
            print(f"{seq:<35} {species:<3} {identity:<6.1f} {species_list}")
    
    # Look for Lactobacillus-specific candidates
    lacto_candidates = clusters_df[
        clusters_df['Species_List'].str.contains('L_', na=False) &
        (clusters_df['Num_Species'] <= 5) &
        (clusters_df['Max_Human_Identity'] < 60)
    ].head(10)
    
    print(f"\nLactobacillus-enriched candidates:")
    if len(lacto_candidates) > 0:
        print(f"{'Sequence':<35} {'Sp':<3} {'ID%':<6} {'Rank':<5} {'Species'}")
        print("-" * 90)
        for _, row in lacto_candidates.iterrows():
            seq = row['Representative_Sequence'][:33]
            species = row['Num_Species']
            identity = row['Max_Human_Identity']
            rank = clusters_df.index.get_loc(row.name) + 1
            species_list = row['Species_List'][:35]
            print(f"{seq:<35} {species:<3} {identity:<6.1f} {rank:<5} {species_list}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("LOW-IDENTITY PRIORITIZATION COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Review TIER1 results")
    print("  2. Check ANATOMIA intensities")
    print("  3. Validate with metaproteome data")
    print("  4. Proceed with alignments of top candidates")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_high_complexity.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    prioritize_clusters_low_identity(input_file, output_dir, log_dir)