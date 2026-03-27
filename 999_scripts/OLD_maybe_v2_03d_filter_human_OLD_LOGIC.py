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
    Uses the SAME logic as the old script for consistency.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_03d_filter_human_OLD_LOGIC")
    sys.stdout = logger
    
    print("="*80)
    print("HUMAN PROTEOME SIMILARITY FILTERING (OLD LOGIC)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    print("\nFiltering Logic (matching old pipeline):")
    print("  Calculate: Alignment_Coverage = (Alignment_Length / Query_Length) * 100")
    print("  EXCLUDE if: identity >= 75% AND coverage >= 70% (both must be true)")
    print("  Then create tiers based on remaining epitopes")
    
    # Load human BLAST results
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    human_blast_file = f"{output_dir}/epitopes_vs_human_ALL_hits.txt"
    print(f"\nLoading: {human_blast_file}")
    
    blast_df = pd.read_csv(human_blast_file, sep='\t')
    print(f"Total BLAST hits: {len(blast_df)}")
    print(f"Columns: {blast_df.columns.tolist()}")
    
    # Calculate alignment coverage (like old script)
    print("\nCalculating Alignment_Coverage...")
    blast_df['Alignment_Coverage'] = (blast_df['Alignment_Length'] / blast_df['Query_Length']) * 100
    
    print(f"Coverage statistics:")
    print(f"  Mean: {blast_df['Alignment_Coverage'].mean():.1f}%")
    print(f"  Median: {blast_df['Alignment_Coverage'].median():.1f}%")
    print(f"  Max: {blast_df['Alignment_Coverage'].max():.1f}%")
    
    # Compare with Query_Coverage from BLAST
    print(f"\nComparing Alignment_Coverage vs Query_Coverage:")
    print(f"  Alignment_Coverage mean: {blast_df['Alignment_Coverage'].mean():.1f}%")
    print(f"  Query_Coverage mean: {blast_df['Query_Coverage'].mean():.1f}%")
    
    # Find epitopes with high human similarity (OLD CRITERIA)
    print("\n" + "="*80)
    print("FILTERING HIGH HUMAN SIMILARITY")
    print("="*80)
    
    high_similarity = blast_df[
        (blast_df['Percent_Identity'] >= 75) & 
        (blast_df['Alignment_Coverage'] >= 70)
    ]
    
    epitopes_to_exclude = set(high_similarity['Epitope_Sequence'].unique())
    
    print(f"\nEpitopes with HIGH human similarity:")
    print(f"  (identity >= 75% AND coverage >= 70%)")
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
    
    # Fill NaN
    epitopes_with_human['Max_Human_Identity'] = epitopes_with_human['Max_Human_Identity'].fillna(0)
    epitopes_with_human['Max_Human_Coverage'] = epitopes_with_human['Max_Human_Coverage'].fillna(0)
    epitopes_with_human['Best_Human_Match'] = epitopes_with_human['Best_Human_Match'].fillna('No_hits')
    
    # Mark bacteria-specific (NOT in exclude list)
    epitopes_with_human['Bacteria_Specific'] = ~epitopes_with_human['Epitope_Sequence'].isin(epitopes_to_exclude)
    
    print(f"\nEpitopes with no human hits: {(epitopes_with_human['Max_Human_Identity'] == 0).sum()}")
    print(f"Epitopes with human hits: {(epitopes_with_human['Max_Human_Identity'] > 0).sum()}")
    print(f"Bacteria-specific (passed filter): {epitopes_with_human['Bacteria_Specific'].sum()}")
    print(f"Excluded (high similarity): {(~epitopes_with_human['Bacteria_Specific']).sum()}")
    
    # Create tiers for bacteria-specific epitopes
    print("\n" + "="*80)
    print("CREATING TIERS")
    print("="*80)
    
    bacteria_specific = epitopes_with_human[epitopes_with_human['Bacteria_Specific'] == True].copy()
    
    def classify_tier(row):
        identity = row['Max_Human_Identity']
        coverage = row['Max_Human_Coverage']
        
        # All of these already passed the exclude filter
        if identity < 60 and coverage < 70:
            return 'BACTERIA_SPECIFIC'
        elif identity < 65 and coverage < 70:
            return 'TIER1'
        elif identity < 70 and coverage < 70:
            return 'TIER2'
        elif identity < 75 and coverage < 70:
            return 'TIER3'
        else:
            # Edge case - shouldn't happen if filter worked
            return 'TIER3'
    
    bacteria_specific['Human_Filter_Category'] = bacteria_specific.apply(classify_tier, axis=1)
    
    # Add category to main dataframe
    epitopes_with_human.loc[epitopes_with_human['Bacteria_Specific'] == False, 'Human_Filter_Category'] = 'EXCLUDE'
    epitopes_with_human.loc[epitopes_with_human['Bacteria_Specific'] == True, 'Human_Filter_Category'] = \
        bacteria_specific['Human_Filter_Category'].values
    
    # Summary
    print("\nClassification Results:")
    category_counts = epitopes_with_human['Human_Filter_Category'].value_counts()
    
    total = len(epitopes_with_human)
    for category in ['BACTERIA_SPECIFIC', 'TIER1', 'TIER2', 'TIER3', 'EXCLUDE']:
        count = category_counts.get(category, 0)
        pct = (count / total) * 100
        print(f"  {category:20} {count:6} epitopes ({pct:5.1f}%)")
    
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
    print("METHOD DISTRIBUTION")
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
    
    logger.close()

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    filter_human_similarity(output_dir, log_dir)