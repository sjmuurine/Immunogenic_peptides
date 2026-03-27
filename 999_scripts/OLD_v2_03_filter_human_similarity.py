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
    logger = Logger(log_dir, "v2_03_human_filter")
    sys.stdout = logger
    
    print("="*80)
    print("HUMAN PROTEOME SIMILARITY FILTERING")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Define thresholds
    thresholds = {
        'BACTERIA_SPECIFIC': {'identity': 60, 'coverage': 70},
        'TIER1': {'identity': 65, 'coverage': 70},
        'TIER2': {'identity': 70, 'coverage': 70},
        'TIER3': {'identity': 75, 'coverage': 70},
        'EXCLUDE': {'identity': 75, 'coverage': 80}
    }
    
    print("\nFiltering Criteria:")
    print(f"  BACTERIA-SPECIFIC: identity < 60% AND coverage < 70%")
    print(f"  TIER1:             identity < 65% AND coverage < 70%")
    print(f"  TIER2:             identity < 70% AND coverage < 70%")
    print(f"  TIER3:             identity < 75% AND coverage < 70%")
    print(f"  EXCLUDE:           identity >= 75% OR coverage >= 80%")
    
    # Load human BLAST results
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    human_blast_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/ARCHIVE_old_species_parsing/003_result_comparisons/epitopes_vs_human_ALL_hits.txt'
    
    print(f"\nLoading human BLAST results from:")
    print(f"  {human_blast_file}")
    
    if not os.path.exists(human_blast_file):
        print(f"\n❌ ERROR: File not found!")
        print(f"Please provide the correct path to human BLAST results.")
        logger.close()
        return
    
    human_df = pd.read_csv(human_blast_file, sep='\t')
    
    print(f"\nLoaded: {len(human_df)} human BLAST hits")
    print(f"Columns: {human_df.columns.tolist()}")
    
    # Check for required columns
    required_cols = ['Epitope_Sequence', 'Percent_Identity', 'Query_Coverage']
    missing = [col for col in required_cols if col not in human_df.columns]
    
    if missing:
        print(f"\n❌ ERROR: Missing required columns: {missing}")
        print(f"Available columns: {human_df.columns.tolist()}")
        logger.close()
        return
    
    # Get best (worst) human match per epitope
    print("\nFinding best human match per epitope...")
    
    best_human = human_df.groupby('Epitope_Sequence').agg({
        'Percent_Identity': 'max',
        'Query_Coverage': 'max',
        'Subject_ID': 'first'  # Keep one example human protein
    }).reset_index()
    
    best_human.columns = ['Epitope_Sequence', 'Max_Human_Identity', 'Max_Human_Coverage', 'Best_Human_Match']
    
    print(f"Unique epitopes with human hits: {len(best_human)}")
    
    # Load all epitopes from comparison (IEDB + k-mer)
    print("\nLoading all epitopes from method comparison...")
    
    comparison_file = f"{output_dir}/IEDB_HIGH_QUALITY_matches_vs_kmer_comparison.txt"
    
    if not os.path.exists(comparison_file):
        print(f"\n❌ ERROR: Comparison file not found: {comparison_file}")
        print("Please run v2_02_compare_iedb_vs_kmer.py first!")
        logger.close()
        return
    
    all_epitopes = pd.read_csv(comparison_file, sep='\t')
    
    print(f"Total epitopes from comparison: {len(all_epitopes)}")
    
    # Merge with human BLAST results
    print("\nMerging epitopes with human similarity data...")
    
    epitopes_with_human = all_epitopes.merge(
        best_human,
        left_on='IEDB_Sequence',
        right_on='Epitope_Sequence',
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
    
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        count = category_counts.get(category, 0)
        pct = (count / len(epitopes_with_human)) * 100
        print(f"  {category:20} {count:6} epitopes ({pct:5.1f}%)")
    
    # Save complete results
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
    
    # Detailed statistics
    print("\n" + "="*80)
    print("DETAILED STATISTICS")
    print("="*80)
    
    print("\nHuman similarity distribution:")
    print(f"  Mean identity: {epitopes_with_human['Max_Human_Identity'].mean():.1f}%")
    print(f"  Median identity: {epitopes_with_human['Max_Human_Identity'].median():.1f}%")
    print(f"  Max identity: {epitopes_with_human['Max_Human_Identity'].max():.1f}%")
    
    print(f"\n  Mean coverage: {epitopes_with_human['Max_Human_Coverage'].mean():.1f}%")
    print(f"  Median coverage: {epitopes_with_human['Max_Human_Coverage'].median():.1f}%")
    print(f"  Max coverage: {epitopes_with_human['Max_Human_Coverage'].max():.1f}%")
    
    # Method distribution per category
    print("\n" + "="*80)
    print("METHOD DISTRIBUTION BY CATEGORY")
    print("="*80)
    
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3']:
        category_df = epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == category]
        
        if len(category_df) > 0:
            print(f"\n{category}:")
            method_counts = category_df['Method_Category'].value_counts()
            for method, count in method_counts.items():
                pct = (count / len(category_df)) * 100
                print(f"  {method:20} {count:5} ({pct:5.1f}%)")
    
    # Summary for publication
    print("\n" + "="*80)
    print("SUMMARY FOR PUBLICATION")
    print("="*80)
    
    total_epitopes = len(epitopes_with_human)
    bacteria_specific = len(epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == 'BACTERIA_SPECIFIC'])
    excluded = len(epitopes_with_human[epitopes_with_human['Human_Filter_Category'] == 'EXCLUDE'])
    
    print(f"\nTotal unique epitopes analyzed: {total_epitopes}")
    print(f"Bacteria-specific (< 60% identity, < 70% coverage): {bacteria_specific} ({bacteria_specific/total_epitopes*100:.1f}%)")
    print(f"Excluded due to human similarity: {excluded} ({excluded/total_epitopes*100:.1f}%)")
    print(f"Retained for further analysis (TIER1-3): {total_epitopes - excluded - bacteria_specific}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("HUMAN FILTERING COMPLETE")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    filter_human_similarity(output_dir, log_dir)