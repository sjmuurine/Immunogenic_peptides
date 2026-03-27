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

def collect_all_epitopes(output_dir, log_dir):
    """
    Collect ALL unique epitope sequences from IEDB, k-mer EXACT, and k-mer FUZZY.
    Track source methods and proteins for each epitope.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_03a_collect_all_epitopes")
    sys.stdout = logger
    
    print("="*80)
    print("COLLECTING ALL EPITOPES FROM ALL METHODS")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    all_epitopes = []
    
    # ========================================================================
    # 1. COLLECT IEDB EPITOPES
    # ========================================================================
    print("\n" + "="*80)
    print("1. LOADING IEDB EPITOPES")
    print("="*80)
    
    iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
    
    print(f"\nLoading: {iedb_file}")
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    
    print(f"Total IEDB BLAST records: {len(iedb_df)}")
    print(f"Unique epitope sequences: {iedb_df['Epitope_Sequence'].nunique()}")
    
    # Group by epitope sequence to get source proteins
    iedb_grouped = iedb_df.groupby('Epitope_Sequence').agg({
        'Protein_Full_Header': lambda x: list(x.unique()),
        'Epitope_ID': 'first',
        'Data_Source': lambda x: '; '.join(x.unique())
    }).reset_index()
    
    for _, row in iedb_grouped.iterrows():
        all_epitopes.append({
            'Epitope_Sequence': row['Epitope_Sequence'],
            'Method': 'IEDB',
            'Source_Proteins': '; '.join(row['Protein_Full_Header'][:10]),  # Limit to 10 for space
            'Num_Source_Proteins': len(row['Protein_Full_Header']),
            'Epitope_ID': row['Epitope_ID'],
            'Data_Source': row['Data_Source']
        })
    
    print(f"✓ Collected {len(iedb_grouped)} unique IEDB epitopes")
    
    # ========================================================================
    # 2. COLLECT K-MER EXACT EPITOPES
    # ========================================================================
    print("\n" + "="*80)
    print("2. LOADING K-MER EXACT EPITOPES")
    print("="*80)
    
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    
    print(f"\nLoading: {exact_file}")
    exact_df = pd.read_csv(exact_file, sep='\t', low_memory=False)
    
    print(f"Total k-mer EXACT motifs: {len(exact_df)}")
    
    for _, row in exact_df.iterrows():
        # Parse protein IDs from Sequence_IDs (comma-separated)
        seq_ids = str(row['Sequence_IDs'])
        protein_list = [pid.strip() for pid in seq_ids.split(',')]
        
        all_epitopes.append({
            'Epitope_Sequence': row['Motif'],
            'Method': 'Kmer_EXACT',
            'Source_Proteins': '; '.join(protein_list[:10]),  # Limit to 10
            'Num_Source_Proteins': row['UniqueSeqCount'],
            'Epitope_ID': f"exact_{row['Motif']}",
            'Data_Source': 'Mixed'
        })
    
    print(f"✓ Collected {len(exact_df)} k-mer EXACT epitopes")
    
    # ========================================================================
    # 3. COLLECT K-MER FUZZY EPITOPES
    # ========================================================================
    print("\n" + "="*80)
    print("3. LOADING K-MER FUZZY EPITOPES")
    print("="*80)
    
    fuzzy_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_fuzzy/epitope_detailed.csv'
    
    print(f"\nLoading: {fuzzy_file}")
    fuzzy_df = pd.read_csv(fuzzy_file)
    
    print(f"Total k-mer FUZZY records: {len(fuzzy_df)}")
    print(f"Unique consensus sequences: {fuzzy_df['Consensus'].nunique()}")
    
    # Group by consensus sequence
    fuzzy_grouped = fuzzy_df.groupby('Consensus').agg({
        'Source_Protein': lambda x: list(x.unique()),
        'Epitope_ID': 'first'
    }).reset_index()
    
    for _, row in fuzzy_grouped.iterrows():
        all_epitopes.append({
            'Epitope_Sequence': row['Consensus'],
            'Method': 'Kmer_FUZZY',
            'Source_Proteins': '; '.join(row['Source_Protein'][:10]),  # Limit to 10
            'Num_Source_Proteins': len(row['Source_Protein']),
            'Epitope_ID': row['Epitope_ID'],
            'Data_Source': 'Mixed'
        })
    
    print(f"✓ Collected {fuzzy_grouped['Consensus'].nunique()} unique k-mer FUZZY epitopes")
    
    # ========================================================================
    # 4. COMBINE AND DEDUPLICATE
    # ========================================================================
    print("\n" + "="*80)
    print("4. COMBINING AND DEDUPLICATING")
    print("="*80)
    
    all_epitopes_df = pd.DataFrame(all_epitopes)
    
    print(f"\nTotal epitope entries (before dedup): {len(all_epitopes_df)}")
    
    # Find duplicates (same sequence in multiple methods)
    duplicates = all_epitopes_df.groupby('Epitope_Sequence').size()
    multi_method = duplicates[duplicates > 1]
    
    print(f"Sequences found in multiple methods: {len(multi_method)}")
    
    if len(multi_method) > 0:
        print("\nTop 10 sequences found in multiple methods:")
        for seq, count in multi_method.head(10).items():
            methods = all_epitopes_df[all_epitopes_df['Epitope_Sequence'] == seq]['Method'].tolist()
            print(f"  {seq}: {methods}")
    
    # Consolidate: merge methods for same sequence
    print("\nConsolidating epitopes...")
    
    consolidated = all_epitopes_df.groupby('Epitope_Sequence').agg({
        'Method': lambda x: '+'.join(sorted(set(x))),
        'Source_Proteins': 'first',  # Take first (could merge if needed)
        'Num_Source_Proteins': 'sum',
        'Epitope_ID': 'first',
        'Data_Source': 'first'
    }).reset_index()
    
    print(f"Unique epitope sequences: {len(consolidated)}")
    
    # Count method combinations
    print("\nMethod distribution:")
    method_counts = consolidated['Method'].value_counts()
    for method, count in method_counts.items():
        pct = (count / len(consolidated)) * 100
        print(f"  {method:30} {count:6} ({pct:5.1f}%)")
    
    # ========================================================================
    # 5. SAVE RESULTS
    # ========================================================================
    print("\n" + "="*80)
    print("5. SAVING RESULTS")
    print("="*80)
    
    # Save complete epitope list
    output_file = f"{output_dir}/all_epitopes_consolidated.txt"
    consolidated.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ Saved consolidated epitopes: {output_file}")
    
    # Save FASTA for BLAST
    fasta_file = f"{output_dir}/all_epitopes_for_blast.fasta"
    
    with open(fasta_file, 'w') as f:
        for idx, row in consolidated.iterrows():
            # Use index as ID for easier tracking
            f.write(f">{idx}\n{row['Epitope_Sequence']}\n")
    
    print(f"✓ Saved FASTA for BLAST: {fasta_file}")
    print(f"  Total sequences: {len(consolidated)}")
    
    # Create ID mapping file
    mapping_file = f"{output_dir}/epitope_id_mapping.txt"
    consolidated[['Epitope_Sequence']].reset_index().to_csv(mapping_file, sep='\t', index=False, header=['ID', 'Epitope_Sequence'])
    print(f"✓ Saved ID mapping: {mapping_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    print(f"\nTotal unique epitope sequences: {len(consolidated)}")
    print(f"Sequences from IEDB only: {(consolidated['Method'] == 'IEDB').sum()}")
    print(f"Sequences from k-mer EXACT only: {(consolidated['Method'] == 'Kmer_EXACT').sum()}")
    print(f"Sequences from k-mer FUZZY only: {(consolidated['Method'] == 'Kmer_FUZZY').sum()}")
    
    # Fix for the backslash issue
    plus_pattern = r'\+'
    multi_method_count = consolidated['Method'].str.contains(plus_pattern).sum()
    print(f"Sequences from multiple methods: {multi_method_count}")
    
    print(f"\nAverage source proteins per epitope: {consolidated['Num_Source_Proteins'].mean():.1f}")
    print(f"Median source proteins per epitope: {consolidated['Num_Source_Proteins'].median():.0f}")
    print(f"Max source proteins per epitope: {consolidated['Num_Source_Proteins'].max()}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("COLLECTION COMPLETE - Ready for BLAST")
    print("="*80)
    
    logger.close()
    
    return fasta_file, output_file

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    
    fasta_file, consolidated_file = collect_all_epitopes(output_dir, log_dir)
    
    print(f"\n✓ FASTA file ready for BLAST: {fasta_file}")
    print(f"✓ Next step: Run BLAST script")