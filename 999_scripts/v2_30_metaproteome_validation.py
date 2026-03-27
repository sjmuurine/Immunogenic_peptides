import pandas as pd
import sys
import os
from datetime import datetime
import subprocess
import tempfile

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

def exact_string_matching(epitope_fasta, metaproteome_fasta, output_dir):
    """
    Perform exact string matching of epitopes in metaproteome sequences.
    """
    print("\n" + "="*60)
    print("EXACT STRING MATCHING")
    print("="*60)
    
    # Load epitope sequences
    epitope_seqs = {}
    with open(epitope_fasta, 'r') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
            elif current_id:
                epitope_seqs[current_id] = line
    
    print(f"Loaded {len(epitope_seqs)} epitope sequences")
    
    # Load metaproteome sequences
    metaproteome_seqs = {}
    with open(metaproteome_fasta, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    metaproteome_seqs[current_id] = ''.join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            metaproteome_seqs[current_id] = ''.join(current_seq)
    
    print(f"Loaded {len(metaproteome_seqs)} metaproteome sequences")
    
    # Perform exact matching
    exact_matches = {}
    
    for epitope_id, epitope_seq in epitope_seqs.items():
        matches = []
        for meta_id, meta_seq in metaproteome_seqs.items():
            if epitope_seq in meta_seq:
                matches.append(meta_id)
        
        exact_matches[epitope_id] = {
            'sequence': epitope_seq,
            'found_exact': len(matches) > 0,
            'num_matches': len(matches),
            'matched_proteins': '; '.join(matches[:10])  # Limit for readability
        }
    
    # Save exact matching results
    exact_results = []
    for epitope_id, data in exact_matches.items():
        exact_results.append({
            'Epitope_ID': epitope_id,
            'Epitope_Sequence': data['sequence'],
            'Found_Exact_Match': data['found_exact'],
            'Num_Exact_Matches': data['num_matches'],
            'Matched_Proteins_Exact': data['matched_proteins']
        })
    
    exact_df = pd.DataFrame(exact_results)
    exact_file = f"{output_dir}/epitopes_exact_matching_results.txt"
    exact_df.to_csv(exact_file, sep='\t', index=False)
    
    found_exact = exact_df['Found_Exact_Match'].sum()
    print(f"\nExact matching results:")
    print(f"  Found exact matches: {found_exact}/{len(exact_df)} ({found_exact/len(exact_df)*100:.1f}%)")
    print(f"  Saved: {exact_file}")
    
    return exact_df

def blast_matching(epitope_fasta, metaproteome_fasta, output_dir, temp_dir):
    """
    Perform BLAST matching of epitopes against metaproteome.
    """
    print("\n" + "="*60)
    print("BLAST MATCHING")
    print("="*60)
    
    # Create BLAST database
    blast_db = f"{temp_dir}/metaproteome_db"
    print(f"Creating BLAST database: {blast_db}")
    
    cmd_makedb = [
        'makeblastdb',
        '-in', metaproteome_fasta,
        '-dbtype', 'prot',
        '-out', blast_db
    ]
    
    try:
        result = subprocess.run(cmd_makedb, capture_output=True, text=True, check=True)
        print("✓ BLAST database created successfully")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating BLAST database: {e}")
        print(f"stderr: {e.stderr}")
        return None
    
    # Run BLAST
    blast_output = f"{temp_dir}/blast_results.txt"
    print(f"Running BLAST search...")
    
    cmd_blast = [
        'blastp',
        '-query', epitope_fasta,
        '-db', blast_db,
        '-out', blast_output,
        '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs',
        '-evalue', '1000',  # Permissive e-value
        '-max_target_seqs', '10',  # Limit hits per query
        '-word_size', '2',  # Sensitive for short sequences
        '-matrix', 'BLOSUM62',
        '-comp_based_stats', '0'  # Turn off composition-based stats for short sequences
    ]
    
    try:
        result = subprocess.run(cmd_blast, capture_output=True, text=True, check=True)
        print("✓ BLAST search completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running BLAST: {e}")
        print(f"stderr: {e.stderr}")
        return None
    
    # Parse BLAST results
    blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs']
    
    try:
        blast_df = pd.read_csv(blast_output, sep='\t', names=blast_columns)
        print(f"BLAST hits found: {len(blast_df)}")
    except:
        print("No BLAST hits found")
        blast_df = pd.DataFrame(columns=blast_columns)
    
    # Summarize BLAST results per epitope
    epitope_blast_results = {}
    
    # Load epitope sequences to get all epitopes (even those without hits)
    epitope_seqs = {}
    with open(epitope_fasta, 'r') as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
            elif current_id:
                epitope_seqs[current_id] = line
    
    # Initialize all epitopes
    for epitope_id, seq in epitope_seqs.items():
        epitope_blast_results[epitope_id] = {
            'sequence': seq,
            'found_blast': False,
            'num_hits': 0,
            'best_identity': 0,
            'best_coverage': 0,
            'best_evalue': float('inf'),
            'matched_proteins': ''
        }
    
    # Process BLAST hits
    if len(blast_df) > 0:
        for epitope_id, group in blast_df.groupby('qseqid'):
            # Get best hit stats
            best_hit = group.loc[group['bitscore'].idxmax()]
            
            epitope_blast_results[epitope_id].update({
                'found_blast': True,
                'num_hits': len(group),
                'best_identity': best_hit['pident'],
                'best_coverage': best_hit['qcovs'] if 'qcovs' in best_hit and pd.notna(best_hit['qcovs']) else 0,
                'best_evalue': best_hit['evalue'],
                'matched_proteins': '; '.join(group['sseqid'].head(5).tolist())
            })
    
    # Create BLAST results dataframe
    blast_results = []
    for epitope_id, data in epitope_blast_results.items():
        blast_results.append({
            'Epitope_ID': epitope_id,
            'Epitope_Sequence': data['sequence'],
            'Found_BLAST_Match': data['found_blast'],
            'Num_BLAST_Hits': data['num_hits'],
            'Best_Identity': data['best_identity'],
            'Best_Coverage': data['best_coverage'],
            'Best_Evalue': data['best_evalue'],
            'Matched_Proteins_BLAST': data['matched_proteins']
        })
    
    blast_results_df = pd.DataFrame(blast_results)
    blast_file = f"{output_dir}/epitopes_blast_matching_results.txt"
    blast_results_df.to_csv(blast_file, sep='\t', index=False)
    
    found_blast = blast_results_df['Found_BLAST_Match'].sum()
    print(f"\nBLAST matching results:")
    print(f"  Found BLAST matches: {found_blast}/{len(blast_results_df)} ({found_blast/len(blast_results_df)*100:.1f}%)")
    if found_blast > 0:
        high_identity = (blast_results_df['Best_Identity'] >= 90).sum()
        print(f"  High identity (≥90%): {high_identity}")
        high_coverage = (blast_results_df['Best_Coverage'] >= 80).sum()
        print(f"  High coverage (≥80%): {high_coverage}")
    print(f"  Saved: {blast_file}")
    
    return blast_results_df

def validate_metaproteome(metaproteome_dir, log_dir):
    """
    Main function to validate epitopes against metaproteome using both exact and BLAST matching.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_30_metaproteome_validation")
    sys.stdout = logger
    
    print("="*80)
    print("METAPROTEOME VALIDATION")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Input files
    epitope_clusters_file = f"{metaproteome_dir}/REFINED_TIER1_clusters.txt"
    epitope_fasta = f"{metaproteome_dir}/refined_tier1_epitopes.fasta"
    metaproteome_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/metaproteome_sequences.fa'
    
    print(f"\nInput files:")
    print(f"  Epitope clusters: {epitope_clusters_file}")
    print(f"  Epitope FASTA: {epitope_fasta}")
    print(f"  Metaproteome FASTA: {metaproteome_fasta}")
    
    # Check files exist
    for file_path in [epitope_clusters_file, epitope_fasta, metaproteome_fasta]:
        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            logger.close()
            return
    
    print("✓ All input files found")
    
    # Load epitope cluster data
    clusters_df = pd.read_csv(epitope_clusters_file, sep='\t')
    print(f"\nLoaded {len(clusters_df)} epitope clusters for validation")
    
    # Create temporary directory for BLAST
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temporary directory: {temp_dir}")
        
        # Perform exact string matching
        exact_results = exact_string_matching(epitope_fasta, metaproteome_fasta, metaproteome_dir)
        
        # Perform BLAST matching
        blast_results = blast_matching(epitope_fasta, metaproteome_fasta, metaproteome_dir, temp_dir)
    
    # Combine results
    if exact_results is not None and blast_results is not None:
        print("\n" + "="*80)
        print("COMBINING RESULTS")
        print("="*80)
        
        # Merge exact and BLAST results
        combined_results = exact_results.merge(blast_results, on=['Epitope_ID', 'Epitope_Sequence'], how='outer')
        
        # Add metaproteome validation to original cluster data
        clusters_with_metaproteome = clusters_df.merge(combined_results, left_on='Cluster_ID', right_on='Epitope_ID', how='left')
        
        # Create summary validation column
        clusters_with_metaproteome['Metaproteome_Detected'] = (
            clusters_with_metaproteome['Found_Exact_Match'].fillna(False) | 
            clusters_with_metaproteome['Found_BLAST_Match'].fillna(False)
        )
        
        clusters_with_metaproteome['Metaproteome_Detection_Type'] = 'Not_Found'
        clusters_with_metaproteome.loc[clusters_with_metaproteome['Found_Exact_Match'] == True, 'Metaproteome_Detection_Type'] = 'Exact_Match'
        clusters_with_metaproteome.loc[
            (clusters_with_metaproteome['Found_Exact_Match'] != True) & 
            (clusters_with_metaproteome['Found_BLAST_Match'] == True), 
            'Metaproteome_Detection_Type'
        ] = 'BLAST_Only'
        
        # Save combined results
        combined_file = f"{metaproteome_dir}/REFINED_TIER1_with_metaproteome_validation.txt"
        clusters_with_metaproteome.to_csv(combined_file, sep='\t', index=False)
        
        # Summary statistics
        total_detected = clusters_with_metaproteome['Metaproteome_Detected'].sum()
        exact_only = (clusters_with_metaproteome['Metaproteome_Detection_Type'] == 'Exact_Match').sum()
        blast_only = (clusters_with_metaproteome['Metaproteome_Detection_Type'] == 'BLAST_Only').sum()
        
        print(f"\nCombined validation results:")
        print(f"  Total epitopes validated: {len(clusters_with_metaproteome)}")
        print(f"  Found in metaproteome: {total_detected}/{len(clusters_with_metaproteome)} ({total_detected/len(clusters_with_metaproteome)*100:.1f}%)")
        print(f"    - Exact matches: {exact_only}")
        print(f"    - BLAST only: {blast_only}")
        print(f"    - Not found: {len(clusters_with_metaproteome) - total_detected}")
        print(f"\n✓ Saved combined results: {combined_file}")
        
        # Show top validated candidates
        validated_candidates = clusters_with_metaproteome[clusters_with_metaproteome['Metaproteome_Detected'] == True]
        if len(validated_candidates) > 0:
            print(f"\nTop 10 metaproteome-validated candidates:")
            print(f"{'Sequence':<35} {'Sp':<3} {'ID%':<6} {'Detection':<12} {'Abundance':<12}")
            print("-" * 85)
            for _, row in validated_candidates.head(10).iterrows():
                seq = row['Representative_Sequence'][:33]
                species = row['Num_Species']
                identity = row['Max_Human_Identity']
                detection = row['Metaproteome_Detection_Type']
                abundance = row.get('Total_Relative_Abundance', 0)
                print(f"{seq:<35} {species:<3} {identity:<6.1f} {detection:<12} {abundance:<12.4f}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("METAPROTEOME VALIDATION COMPLETE")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    metaproteome_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    validate_metaproteome(metaproteome_dir, log_dir)