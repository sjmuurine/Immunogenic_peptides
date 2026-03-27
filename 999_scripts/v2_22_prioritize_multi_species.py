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
    Calculate priority score emphasizing low human identity AND multi-species presence.
    
    Score = (Species × 50) + (IEDB × 30) - (Identity × 20.0) + (Low_Identity_Bonus) - (Single_Species_Penalty)
    
    Higher is better.
    Strong emphasis on bacteria-specific epitopes with broad species coverage.
    """
    species_score = row['Num_Species'] * 50  # Base species score
    iedb_bonus = 30 if row['Has_IEDB'] else 0
    identity_penalty = row['Max_Human_Identity'] * 20.0
    
    # Extra bonus for very low identity (bacteria-specific)
    identity = row['Max_Human_Identity']
    low_identity_bonus = 0
    if identity < 30:
        low_identity_bonus = 200
    elif identity < 50:
        low_identity_bonus = 100
    elif identity < 60:
        low_identity_bonus = 50
    
    # Penalty for single-species epitopes (less useful for broad detection)
    single_species_penalty = 150 if row['Num_Species'] == 1 else 0
    
    total = species_score + iedb_bonus - identity_penalty + low_identity_bonus - single_species_penalty
    
    return total

def prioritize_clusters_multi_species(input_file, output_dir, log_dir):
    """
    Prioritize clusters with emphasis on low human identity and multi-species presence.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_22_prioritize_multi_species")
    sys.stdout = logger
    
    print("="*80)
    print("EPITOPE PRIORITIZATION (LOW IDENTITY + MULTI-SPECIES)")
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
    
    print("\nPriority formula (MULTI-SPECIES FOCUS):")
    print("  Score = (Species × 50) + (IEDB × 30) - (Identity × 20.0) + Low_Identity_Bonus - Single_Species_Penalty")
    print("  Where:")
    print("    Low_Identity_Bonus:")
    print("      - <30% identity: +200 (very bacteria-specific)")
    print("      - 30-50% identity: +100 (bacteria-specific)")
    print("      - 50-60% identity: +50 (moderately specific)")
    print("    Single_Species_Penalty: -150 (pushes 1-species epitopes to TIER2)")
    
    clusters_df['Priority_Score'] = clusters_df.apply(calculate_priority_score, axis=1)
    
    # Sort by priority
    clusters_df = clusters_df.sort_values('Priority_Score', ascending=False)
    
    print(f"\nPriority score statistics:")
    print(f"  Mean: {clusters_df['Priority_Score'].mean():.1f}")
    print(f"  Median: {clusters_df['Priority_Score'].median():.1f}")
    print(f"  Min: {clusters_df['Priority_Score'].min():.1f}")
    print(f"  Max: {clusters_df['Priority_Score'].max():.1f}")
    
    # Show species distribution impact
    print(f"\nSingle-species penalty impact:")
    single_species = (clusters_df['Num_Species'] == 1).sum()
    print(f"  Single-species epitopes: {single_species:,} (penalty: -150)")
    print(f"  Multi-species epitopes: {len(clusters_df) - single_species:,} (no penalty)")
    
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
    
    # Check if single-species epitopes were successfully pushed out
    print("\n" + "="*80)
    print("SINGLE-SPECIES ANALYSIS")
    print("="*80)
    
    single_in_top50 = (top50['Num_Species'] == 1).sum()
    print(f"\nSingle-species epitopes in TOP 50: {single_in_top50}")
    
    if single_in_top50 > 0:
        print("WARNING: Some single-species epitopes still in TOP 50!")
        single_examples = top50[top50['Num_Species'] == 1]
        for _, row in single_examples.iterrows():
            rank = clusters_df.index.get_loc(row.name) + 1
            print(f"  Rank {rank}: {row['Representative_Sequence'][:30]} - {row['Species_List']} (Score: {row['Priority_Score']:.1f})")
    else:
        print("✓ Successfully pushed single-species epitopes out of TOP 50!")
    
    print(f"\nTOP 50 characteristics:")
    print(f"  Species range: {top50['Num_Species'].min()}-{top50['Num_Species'].max()}")
    print(f"  Mean species: {top50['Num_Species'].mean():.1f}")
    print(f"  Identity range: {top50['Max_Human_Identity'].min():.1f}-{top50['Max_Human_Identity'].max():.1f}%")
    print(f"  Mean identity: {top50['Max_Human_Identity'].mean():.1f}%")
    print(f"  IEDB sequences: {top50['Has_IEDB'].sum()}")
    
    print(f"\nSpecies distribution in TOP 50:")
    species_dist = top50['Num_Species'].value_counts().sort_index(ascending=False)
    for species_count, freq in species_dist.items():
        print(f"  {species_count} species: {freq} epitopes")
    
    # Create simplified tiering (2 tiers)
    print("\n" + "="*80)
    print("CREATING TIERS")
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
    
    # Check single-species distribution across tiers
    print(f"\nSingle-species distribution by tier:")
    for tier in ['TIER1', 'TIER2', 'TIER3']:
        tier_df = clusters_df[clusters_df['Tier'] == tier]
        if len(tier_df) > 0:
            single_count = (tier_df['Num_Species'] == 1).sum()
            pct = single_count / len(tier_df) * 100 if len(tier_df) > 0 else 0
            print(f"  {tier}: {single_count}/{len(tier_df)} ({pct:.1f}%) single-species")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    # Save all prioritized clusters
    all_output = f"{output_dir}/epitope_clusters_prioritized_multi_species.txt"
    clusters_df.to_csv(all_output, sep='\t', index=False)
    print(f"\n✓ Saved all prioritized: {all_output}")
    
    # Save by tier
    for tier in ['TIER1', 'TIER2']:
        tier_df = clusters_df[clusters_df['Tier'] == tier]
        if len(tier_df) > 0:
            tier_file = f"{output_dir}/{tier}_clusters_multi_species.txt"
            tier_df.to_csv(tier_file, sep='\t', index=False)
            print(f"✓ Saved {tier}: {tier_file} ({len(tier_df)} clusters)")
    
    # Show examples of pushed-down single-species epitopes
    print("\n" + "="*80)
    print("SINGLE-SPECIES EPITOPES (pushed to TIER2/3)")
    print("="*80)
    
    single_species_all = clusters_df[clusters_df['Num_Species'] == 1]
    if len(single_species_all) > 0:
        print(f"\nTop 5 single-species epitopes (now in lower tiers):")
        single_top5 = single_species_all.head(5)
        for _, row in single_top5.iterrows():
            rank = clusters_df.index.get_loc(row.name) + 1
            seq = row['Representative_Sequence'][:30]
            species = row['Species_List']
            identity = row['Max_Human_Identity']
            tier = row['Tier']
            score = row['Priority_Score']
            print(f"  Rank {rank:3} ({tier}): {seq:30} - {species} ({identity:.1f}% ID, Score: {score:.1f})")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("MULTI-SPECIES PRIORITIZATION COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Review TIER1 results (should have multi-species epitopes)")
    print("  2. Check ANATOMIA intensities")
    print("  3. Validate with metaproteome data")
    print("  4. Proceed with alignments of top candidates")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_high_complexity.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    prioritize_clusters_multi_species(input_file, output_dir, log_dir)