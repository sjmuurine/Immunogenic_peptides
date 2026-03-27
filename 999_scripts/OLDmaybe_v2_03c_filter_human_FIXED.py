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
    
    Criteria:
    - BACTERIA-SPECIFIC: identity < 60% AND coverage < 70%
    - TIER1: identity < 65% AND coverage < 70%
    - TIER2: identity < 70% AND coverage < 70%
    - TIER3: identity < 75% AND coverage < 70%
    - EXCLUDE: identity >= 75% OR coverage >= 80%
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_03c_filter_human_FIXED")
    sys.stdout = logger
    
    print("="*80)
    print("HUMAN PROTEOME SIMILARITY FILTERING (FIXED)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nFiltering Criteria:")
    print(f"  BACTERIA-SPECIFIC: identity < 60% AND coverage < 70%")
    print(f"  TIER1:             identity < 65% AND coverage < 70%")
    print(f"  TIER2:             identity < 70% AND coverage < 70%")
    print(f"  TIER3:             identity < 75% AND coverage < 70%")
    print(f"  EXCLUDE:           identity >= 75% OR coverage >= 80%")
    
    # Load data
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    # Load human BLAST results - FIXED PATH
    human_blast_file = f"{output_dir}/epitopes_vs_human_ALL_hits.txt"
    print(f"\nLoading human BLAST results: {human_blast_file}")
    
    if not os.path.exists(human_blast_file):
        print(f"\n❌ ERROR: File not found: {human_blast_file}")
        logger.close()
        return
    
    human_df = pd.read_csv(human_blast_file, sep='\t')
    print(f"Loaded: {len(human_df)} human BLAST hits")
    print(f"Columns: {human_df.columns.tolist()}")
    
    # Get best (worst) human match per epitope
    print("\nFinding best (worst) human match per epitope...")
    print("  Taking MAX identity and MAX coverage per epitope")
    
    best_human = human_df.groupby('Epitope_Sequence').agg({
        'Percent_Identity': 'max',
        'Query_Coverage': 'max',
        'Subject_ID': 'first'
    }).reset_index()
    
    best_human.columns = ['Epitope_Sequence', 'Max_Human_Identity', 'Max_Human_Coverage', 'Best_Human_Match']
    
    print(f"Unique epitopes with human hits: {len(best_human)}")
    
    # Show some statistics on human matches
    print(f"\nHuman similarity statistics:")
    print(f"  Mean identity: {best_human['Max_Human_Identity'].mean():.1f}%")
    print(f"  Median identity: {best_human['Max_Human_Identity'].median():.1f}%")
    print(f"  Max identity: {best_human['Max_Human_Identity'].max():.1f}%")
    print(f"  Mean coverage: {best_human['Max_Human_Coverage'].mean():.1f}%")
    print(f"  Median coverage: {best_human['Max_Human_Coverage'].median():.1f}%")
    print(f"  Max coverage: {best_human['Max_Human_Coverage'].max():.1f}%")
    
    # Load all epitopes
    all_epitopes_file = f"{output_dir}/all_epitopes_consolidated.txt"
    print(f"\nLoading all epitopes: {all_epitopes_file}")
    
    all_epitopes = pd.read_csv(all_epitopes_file, sep='\t')
    print(f"Total epitopes: {len(all_epitopes)}")
    
    # Merge with human BLAST results
    print("\nMerging epitopes with human similarity data...")
    
    epitopes_with_human = all_epitopes.merge(
        best_human,
        on='Epitope_Sequence',
        how='left'
    )
    
    # Fill NaN (no human hit) with 0
    epitopes_with_human['Max_Human_Identity'] = epitopes_with_human['Max_Human_Identity'].fillna(0)
    epitopes_with_human['Max_Human_Coverage'] = epitopes_with_human['Max_Human_Coverage'].fillna(0)
    epitopes_with_human['Best_Human_Match'] = epitopes_with_human['Best_Human_Match'].fillna('No_hits')
    
    print(f"Epitopes with no human hits: {(epitopes_with_human['Max_Human_Identity'] == 0).sum()}")
    print(f"Epitopes with human hits: {(epitopes_with_human['Max_Human_Identity'] > 0).sum()}")
    
    # Classify epitopes
    print("\n" + "="*80)
    print("CLASSIFYING EPITOPES")
    print("="*80)
    
    def classify_epitope(row):
        identity = row['Max_Human_Identity']
        coverage = row['Max_Human_Coverage']
        
        # EXCLUDE first (most stringent)
        if identity >= 75 or coverage >= 80:
            return 'EXCLUDE'
        
        # Then check tiers (from strictest to most permissive)
        if identity < 60 and coverage < 70:
            return 'BACTERIA_SPECIFIC'
        elif identity < 65 and coverage < 70:
            return 'TIER1'
        elif identity < 70 and coverage < 70:
            return 'TIER2'
        elif identity < 75 and coverage < 70:
            return 'TIER3'
        else:
            return 'EXCLUDE'
    
    epitopes_with_human['Human_Filter_Category'] = epitopes_with_human.apply(classify_epitope, axis=1)
    
    # Summary statistics
    print("\nClassification Results:")
    category_counts = epitopes_with_human['Human_Filter_Category'].value_counts()
    
    total = len(epitopes_with_human)
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        count = category_counts.get(category, 0)
        pct = (count / total) * 100
        print(f"  {category:20} {count:6} epitopes ({pct:5.1f}%)")
    
    # Detailed breakdown of EXCLUDE reasons
    print("\nEXCLUDE breakdown:")
    excluded = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == 'EXCLUDE']
    high_identity = (excluded['Max_Human_Identity'] >= 75).sum()
    high_coverage = (excluded['Max_Human_Coverage'] >= 80).sum()
    both = ((excluded['Max_Human_Identity'] >= 75) & (excluded['Max_Human_Coverage'] >= 80)).sum()
    print(f"  High identity (>=75%): {high_identity}")
    print(f"  High coverage (>=80%): {high_coverage}")
    print(f"  Both: {both}")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    all_results_file = f"{output_dir}/all_epitopes_with_human_similarity.txt"
    epitopes_with_human.to_csv(all_results_file, sep='\t', index=False)
    print(f"\n✓ Saved all epitopes: {all_results_file}")
    
    # Save category-specific files
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        category_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        
        if len(category_df) > 0:
            category_file = f"{output_dir}/{category}_epitopes.txt"
            category_df.to_csv(category_file, sep='\t', index=False)
            print(f"✓ Saved {category}: {category_file} ({len(category_df)} epitopes)")
    
    # Detailed statistics by method
    print("\n" + "="*80)
    print("METHOD DISTRIBUTION BY CATEGORY")
    print("="*80)
    
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3']:
        category_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        
        if len(category_df) > 0:
            print(f"\n{category}:")
            method_counts = category_df['Method'].value_counts()
            for method, count in method_counts.items():
                pct = (count / len(category_df)) * 100
                print(f"  {method:20} {count:5} ({pct:5.1f}%)")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("FILTERING COMPLETE")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    filter_human_similarity(output_dir, log_dir)