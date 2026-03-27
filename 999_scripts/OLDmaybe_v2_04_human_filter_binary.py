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

def filter_human_similarity_binary(output_dir, log_dir):
    """
    Binary human filtering: PASS or EXCLUDE only.
    Tiering will be done later based on multiple factors.
    
    EXCLUDE if: identity >= 75% AND coverage >= 70%
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_04_human_filter_binary")
    sys.stdout = logger
    
    print("="*80)
    print("HUMAN PROTEOME SIMILARITY FILTER (BINARY)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nFiltering Criteria:")
    print("  EXCLUDE: identity >= 75% AND coverage >= 70%")
    print("  PASS: All others")
    print("\nNote: Tiering will be done later based on species, method, etc.")
    
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
    
    # Find epitopes to EXCLUDE
    print("\n" + "="*80)
    print("IDENTIFYING HIGH HUMAN SIMILARITY")
    print("="*80)
    
    high_similarity = blast_df[
        (blast_df['Percent_Identity'] >= 75) & 
        (blast_df['Alignment_Coverage'] >= 70)
    ]
    
    epitopes_to_exclude = set(high_similarity['Epitope_Sequence'].unique())
    
    print(f"\nEpitopes to EXCLUDE (high human similarity):")
    print(f"  identity >= 75% AND coverage >= 70%")
    print(f"  Total: {len(epitopes_to_exclude)}")
    
    # Get best human match for ALL epitopes
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
    
    # Mark PASS vs EXCLUDE
    epitopes_with_human['Passes_Human_Filter'] = ~epitopes_with_human['Epitope_Sequence'].isin(epitopes_to_exclude)
    
    # Summary
    print("\n" + "="*80)
    print("FILTERING RESULTS")
    print("="*80)
    
    passed = epitopes_with_human['Passes_Human_Filter'].sum()
    excluded = (~epitopes_with_human['Passes_Human_Filter']).sum()
    total = len(epitopes_with_human)
    
    print(f"\nTotal epitopes: {total}")
    print(f"  PASS (bacteria-specific): {passed} ({passed/total*100:.1f}%)")
    print(f"  EXCLUDE (too similar): {excluded} ({excluded/total*100:.1f}%)")
    
    # Statistics on passed epitopes
    print("\n" + "="*80)
    print("PASSED EPITOPES - STATISTICS")
    print("="*80)
    
    passed_df = epitopes_with_human[epitopes_with_human['Passes_Human_Filter'] == True]
    
    print(f"\nMethod distribution:")
    method_counts = passed_df['Method'].value_counts()
    for method, count in method_counts.items():
        pct = (count / len(passed_df)) * 100
        print(f"  {method:30} {count:6} ({pct:5.1f}%)")
    
    # Identity/coverage stats for passed
    passed_with_hits = passed_df[passed_df['Max_Human_Identity'] > 0]
    
    if len(passed_with_hits) > 0:
        print(f"\nHuman similarity stats (passed epitopes with hits):")
        print(f"  Mean identity: {passed_with_hits['Max_Human_Identity'].mean():.1f}%")
        print(f"  Median identity: {passed_with_hits['Max_Human_Identity'].median():.1f}%")
        print(f"  Max identity: {passed_with_hits['Max_Human_Identity'].max():.1f}%")
        print(f"  Mean coverage: {passed_with_hits['Max_Human_Coverage'].mean():.1f}%")
        print(f"  Median coverage: {passed_with_hits['Max_Human_Coverage'].median():.1f}%")
        print(f"  Max coverage: {passed_with_hits['Max_Human_Coverage'].max():.1f}%")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    all_results_file = f"{output_dir}/all_epitopes_human_filtered.txt"
    epitopes_with_human.to_csv(all_results_file, sep='\t', index=False)
    print(f"\n✓ Saved all epitopes: {all_results_file}")
    
    # Save passed epitopes separately
    passed_file = f"{output_dir}/epitopes_PASSED_human_filter.txt"
    passed_df.to_csv(passed_file, sep='\t', index=False)
    print(f"✓ Saved PASSED: {passed_file} ({len(passed_df)} epitopes)")
    
    # Save excluded epitopes
    excluded_df = epitopes_with_human[epitopes_with_human['Passes_Human_Filter'] == False]
    excluded_file = f"{output_dir}/epitopes_EXCLUDED_human_similar.txt"
    excluded_df.to_csv(excluded_file, sep='\t', index=False)
    print(f"✓ Saved EXCLUDED: {excluded_file} ({len(excluded_df)} epitopes)")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("BINARY FILTERING COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Add corrected species counts")
    print("  2. Cluster similar epitopes")
    print("  3. Calculate priority scores")
    print("  4. Create tiers based on composite score")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    filter_human_similarity_binary(output_dir, log_dir)