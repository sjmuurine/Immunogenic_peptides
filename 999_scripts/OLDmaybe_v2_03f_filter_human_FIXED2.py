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

def filter_human_similarity(output_dir, log_dir):
    """
    Filter epitopes based on human proteome similarity.
    
    UPDATED CRITERIA:
    - EXCLUDE: identity >= 75% OR coverage >= 80%
    - For remaining epitopes, tier by IDENTITY only:
      - BACTERIA_SPECIFIC: identity < 60%
      - TIER1: identity >= 60% and < 65%
      - TIER2: identity >= 65% and < 70%
      - TIER3: identity >= 70% and < 75%
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_03f_filter_human_FIXED2")
    sys.stdout = logger
    
    print("="*80)
    print("HUMAN PROTEOME SIMILARITY FILTERING (FIXED)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nFiltering Criteria:")
    print("  Step 1 - EXCLUDE: identity >= 75% OR coverage >= 80%")
    print("  Step 2 - Tier remaining by IDENTITY:")
    print("    BACTERIA_SPECIFIC: identity < 60%")
    print("    TIER1:             60% <= identity < 65%")
    print("    TIER2:             65% <= identity < 70%")
    print("    TIER3:             70% <= identity < 75%")
    
    # Load human BLAST results
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    human_blast_file = f"{output_dir}/epitopes_vs_human_ALL_hits.txt"
    print(f"\nLoading: {human_blast_file}")
    
    blast_df = pd.read_csv(human_blast_file, sep='\t')
    print(f"Total BLAST hits: {len(blast_df)}")
    
    # Calculate alignment coverage
    print("\nCalculating Alignment_Coverage...")
    blast_df['Alignment_Coverage'] = (blast_df['Alignment_Length'] / blast_df['Query_Length']) * 100
    
    # Get best human match for ALL epitopes FIRST
    print("\nFinding best human match per epitope...")
    
    best_human = blast_df.groupby('Epitope_Sequence').agg({
        'Percent_Identity': 'max',
        'Alignment_Coverage': 'max',
        'Subject_ID': 'first'
    }).reset_index()
    
    best_human.columns = ['Epitope_Sequence', 'Max_Human_Identity', 'Max_Human_Coverage', 'Best_Human_Match']
    
    print(f"Unique epitopes with human hits: {len(best_human)}")
    
    # Load all epitopes
    all_epitopes_file = f"{output_dir}/all_epitopes_consolidated.txt"
    print(f"\nLoading all epitopes: {all_epitopes_file}")
    
    all_epitopes = pd.read_csv(all_epitopes_file, sep='\t')
    print(f"Total epitopes: {len(all_epitopes)}")
    
    # Merge with human data
    epitopes_with_human = all_epitopes.merge(
        best_human,
        on='Epitope_Sequence',
        how='left'
    )
    
    # Fill NaN (no human hits)
    epitopes_with_human['Max_Human_Identity'] = epitopes_with_human['Max_Human_Identity'].fillna(0)
    epitopes_with_human['Max_Human_Coverage'] = epitopes_with_human['Max_Human_Coverage'].fillna(0)
    epitopes_with_human['Best_Human_Match'] = epitopes_with_human['Best_Human_Match'].fillna('No_hits')
    
    print(f"\nEpitopes with no human hits: {(epitopes_with_human['Max_Human_Identity'] == 0).sum()}")
    print(f"Epitopes with human hits: {(epitopes_with_human['Max_Human_Identity'] > 0).sum()}")
    
    # Classify epitopes
    print("\n" + "="*80)
    print("CLASSIFYING EPITOPES")
    print("="*80)
    
    def classify_epitope(row):
        identity = row['Max_Human_Identity']
        coverage = row['Max_Human_Coverage']
        
        # Step 1: Check if should be EXCLUDED (high similarity)
        if identity >= 75 or coverage >= 80:
            return 'EXCLUDE'
        
        # Step 2: Tier remaining by identity
        if identity < 60:
            return 'BACTERIA_SPECIFIC'
        elif identity < 65:
            return 'TIER1'
        elif identity < 70:
            return 'TIER2'
        else:  # identity < 75 (since >= 75 already excluded)
            return 'TIER3'
    
    epitopes_with_human['Human_Filter_Category'] = epitopes_with_human.apply(classify_epitope, axis=1)
    
    # Summary statistics
    print("\nClassification Results:")
    category_counts = epitopes_with_human['Human_Filter_Category'].value_counts()
    
    total = len(epitopes_with_human)
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        count = category_counts.get(category, 0)
        pct = (count / total) * 100
        print(f"  {category:20} {count:6} epitopes ({pct:5.1f}%)")
    
    # Exclusion breakdown
    print("\n" + "="*80)
    print("EXCLUSION BREAKDOWN")
    print("="*80)
    
    excluded = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == 'EXCLUDE']
    
    high_id_only = ((excluded['Max_Human_Identity'] >= 75) & (excluded['Max_Human_Coverage'] < 80)).sum()
    high_cov_only = ((excluded['Max_Human_Identity'] < 75) & (excluded['Max_Human_Coverage'] >= 80)).sum()
    both = ((excluded['Max_Human_Identity'] >= 75) & (excluded['Max_Human_Coverage'] >= 80)).sum()
    
    print(f"\nExcluded epitopes: {len(excluded)}")
    print(f"  High identity only (>=75%): {high_id_only}")
    print(f"  High coverage only (>=80%): {high_cov_only}")
    print(f"  Both conditions: {both}")
    
    # Identity distribution by category
    print("\n" + "="*80)
    print("IDENTITY DISTRIBUTION BY CATEGORY")
    print("="*80)
    
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3']:
        cat_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        if len(cat_df) > 0:
            ids = cat_df[cat_df['Max_Human_Identity'] > 0]['Max_Human_Identity']
            if len(ids) > 0:
                print(f"\n{category}: {len(cat_df)} epitopes")
                print(f"  Mean identity: {ids.mean():.1f}%")
                print(f"  Median identity: {ids.median():.1f}%")
                print(f"  Min identity: {ids.min():.1f}%")
                print(f"  Max identity: {ids.max():.1f}%")
            else:
                print(f"\n{category}: {len(cat_df)} epitopes (all have 0% identity)")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    all_results_file = f"{output_dir}/all_epitopes_with_human_similarity.txt"
    epitopes_with_human.to_csv(all_results_file, sep='\t', index=False)
    print(f"\n✓ Saved: {all_results_file}")
    
    # Save by category
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        category_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        if len(category_df) > 0:
            category_file = f"{output_dir}/{category}_epitopes.txt"
            category_df.to_csv(category_file, sep='\t', index=False)
            print(f"✓ Saved {category}: {category_file} ({len(category_df)} epitopes)")
    
    # Method distribution
    print("\n" + "="*80)
    print("METHOD DISTRIBUTION BY CATEGORY")
    print("="*80)
    
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3']:
        category_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        if len(category_df) > 0:
            print(f"\n{category}: {len(category_df)} epitopes")
            method_counts = category_df['Method'].value_counts()
            for method, count in method_counts.items():
                pct = (count / len(category_df)) * 100
                print(f"  {method:30} {count:6} ({pct:5.1f}%)")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("FILTERING COMPLETE")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    filter_human_similarity(output_dir, log_dir)