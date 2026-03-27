import pandas as pd
import numpy as np
from collections import defaultdict
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

def parse_species_from_protein_id(protein_id):
    """Parse species from protein ID - all 17 species."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    species_patterns = {
        'fragilis': 'B_fragilis', 'nucleatum': 'F_nucleatum', 'muciniphila': 'A_muciniphila',
        'prausnitzii': 'F_prausnitzii', 'reuteri': 'L_reuteri', 'cloacae': 'E_cloacae',
        'faecalis': 'E_faecalis', 'salivarius': 'S_salivarius', 'longum': 'B_longum',
        'plantarum': 'L_plantarum', 'casei': 'L_casei', 'sanguinis': 'T_sanguinis',
        'magna': 'V_magna', 'russellii': 'P_russellii', 'freudenreichii': 'P_freudenreichii',
        'fermentum': 'L_fermentum', 'adiacens': 'G_adiacens', 'gasseri': 'L_gasseri'
    }
    
    if 'LITERATURE' in protein_id:
        for pattern, species in species_patterns.items():
            if pattern.lower() in protein_id.lower():
                return species
        return None
    
    if ' ' not in protein_id:
        return None
    
    parts = protein_id.split(' ', 1)[1].split('_')
    og_idx = None
    for i, part in enumerate(parts):
        if part.startswith('OG') or part == 'unique':
            og_idx = i
            break
    
    if og_idx and og_idx > 0:
        return '_'.join(parts[:og_idx])
    
    return None

def calculate_simple_overlap_similarity(seq1, seq2):
    """Calculate similarity by trying all possible alignments."""
    if seq1 == seq2:
        return 100.0
    
    max_similarity = 0
    
    for offset in range(-len(seq2) + 1, len(seq1)):
        matches = 0
        overlap_length = 0
        
        for i in range(len(seq1)):
            j = i - offset
            if 0 <= j < len(seq2):
                overlap_length += 1
                if seq1[i] == seq2[j]:
                    matches += 1
        
        if overlap_length > 0:
            min_len = min(len(seq1), len(seq2))
            if overlap_length >= min_len * 0.8:
                similarity = (matches / overlap_length) * 100
                max_similarity = max(max_similarity, similarity)
    
    return max_similarity

def is_substring_match(seq1, seq2, threshold=0.8):
    """Quick check if one sequence is substring of another."""
    if seq1 in seq2 or seq2 in seq1:
        return True
    
    # Check if substantial overlap exists (quick pre-filter)
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    # If length difference > 20%, probably not similar enough
    if abs(len(seq1) - len(seq2)) > max_len * 0.2:
        return False
    
    return None  # Need full check

def cluster_sequences_fast(sequences_data, similarity_threshold=80):
    """Fast clustering using length-based pre-filtering and substring detection."""
    
    n = len(sequences_data)
    print(f"\n  Fast clustering {n} sequences with {similarity_threshold}% threshold...")
    
    # Group by length
    print("  Step 1: Grouping by length...")
    length_groups = defaultdict(list)
    for i, data in enumerate(sequences_data):
        length = len(data['Epitope_Sequence'])
        length_groups[length].append(i)
    
    print(f"    Found {len(length_groups)} unique lengths (range: {min(length_groups.keys())}-{max(length_groups.keys())} aa)")
    
    # Build comparison pairs (only nearby lengths)
    print("  Step 2: Building comparison pairs...")
    comparison_pairs = []
    
    for length in sorted(length_groups.keys()):
        indices = length_groups[length]
        
        # Compare within same length
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                comparison_pairs.append((indices[i], indices[j]))
        
        # Compare with nearby lengths (±1, ±2)
        for offset in [1, 2]:
            if length + offset in length_groups:
                for i in indices:
                    for j in length_groups[length + offset]:
                        if i < j:
                            comparison_pairs.append((i, j))
    
    print(f"    Total comparison pairs: {len(comparison_pairs):,} (vs {n*(n-1)//2:,} without optimization)")
    print(f"    Reduction: {(1 - len(comparison_pairs)/(n*(n-1)//2))*100:.1f}%")
    
    # Perform clustering
    print("  Step 3: Clustering...")
    clusters = []
    assigned = set()
    comparisons_done = 0
    substring_matches = 0
    full_calcs = 0
    
    # Sort indices to process in order
    all_indices = list(range(n))
    
    for i in all_indices:
        if i in assigned:
            continue
        
        if len(clusters) % 5000 == 0 and len(clusters) > 0:
            print(f"    Progress: {len(clusters):,} clusters, {comparisons_done:,} comparisons, {len(assigned):,}/{n} assigned", end='\r')
        
        # Start new cluster
        cluster = [i]
        assigned.add(i)
        seq_i = sequences_data[i]['Epitope_Sequence']
        
        # Find similar sequences
        for idx_i, idx_j in comparison_pairs:
            if idx_i != i:
                continue
            
            j = idx_j
            
            if j in assigned:
                continue
            
            seq_j = sequences_data[j]['Epitope_Sequence']
            
            # Quick substring check
            substring_result = is_substring_match(seq_i, seq_j)
            
            if substring_result == True:
                # Definitely similar
                cluster.append(j)
                assigned.add(j)
                substring_matches += 1
                comparisons_done += 1
            elif substring_result == False:
                # Definitely not similar enough
                comparisons_done += 1
            else:
                # Need full calculation
                similarity = calculate_simple_overlap_similarity(seq_i, seq_j)
                full_calcs += 1
                comparisons_done += 1
                
                if similarity >= similarity_threshold:
                    cluster.append(j)
                    assigned.add(j)
        
        clusters.append(cluster)
    
    print()
    print(f"\n  Clustering stats:")
    print(f"    Total comparisons: {comparisons_done:,}")
    print(f"    Substring matches: {substring_matches:,}")
    print(f"    Full similarity calcs: {full_calcs:,}")
    print(f"    Clusters created: {len(clusters):,}")
    
    return clusters

def cluster_epitopes_fast(input_file, output_dir, log_dir, similarity_threshold=80):
    """
    Fast epitope clustering with species aggregation.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_07_cluster_fast")
    sys.stdout = logger
    
    print("="*80)
    print("FAST EPITOPE CLUSTERING WITH SPECIES AGGREGATION")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Similarity threshold: {similarity_threshold}%")
    
    # Load epitopes
    print("\n" + "="*80)
    print("LOADING DATA")
    print("="*80)
    
    print(f"\nLoading: {input_file}")
    epitopes_df = pd.read_csv(input_file, sep='\t')
    
    print(f"Total epitopes: {len(epitopes_df)}")
    
    # Filter by length
    print("\nFiltering by length (11-24 aa)...")
    epitopes_df = epitopes_df[
        (epitopes_df['Epitope_Sequence'].str.len() >= 11) & 
        (epitopes_df['Epitope_Sequence'].str.len() <= 24)
    ]
    
    print(f"After length filter: {len(epitopes_df)}")
    
    # Prepare data
    sequences_data = epitopes_df.to_dict('records')
    
    # Perform clustering
    print("\n" + "="*80)
    print("CLUSTERING")
    print("="*80)
    
    clusters = cluster_sequences_fast(sequences_data, similarity_threshold)
    
    singleton_clusters = [c for c in clusters if len(c) == 1]
    multi_member_clusters = [c for c in clusters if len(c) > 1]
    
    print(f"\n  Total clusters: {len(clusters):,}")
    print(f"  Singleton clusters: {len(singleton_clusters):,}")
    print(f"  Multi-member clusters: {len(multi_member_clusters):,}")
    
    if multi_member_clusters:
        sizes = [len(c) for c in multi_member_clusters]
        print(f"  Largest cluster: {max(sizes)} members")
        print(f"  Average multi-cluster size: {np.mean(sizes):.1f} members")
    
    # Aggregate metadata
    print("\n" + "="*80)
    print("AGGREGATING METADATA AND PARSING SPECIES")
    print("="*80)
    
    cluster_results = []
    detailed_clusters = []
    
    for cluster_id, cluster_indices in enumerate(clusters, 1):
        if cluster_id % 10000 == 0:
            print(f"  Processing cluster {cluster_id:,}/{len(clusters):,}...", end='\r')
        
        cluster_epitopes = [sequences_data[i] for i in cluster_indices]
        
        # Choose representative (IEDB priority)
        best_idx = cluster_indices[0]
        best_method = sequences_data[best_idx].get('Method', '')
        best_has_iedb = 'IEDB' in best_method
        
        for idx in cluster_indices[1:]:
            method = sequences_data[idx].get('Method', '')
            has_iedb = 'IEDB' in method
            
            if has_iedb and not best_has_iedb:
                best_idx = idx
                best_method = method
                best_has_iedb = True
            elif has_iedb == best_has_iedb:
                if sequences_data[idx].get('Num_Source_Proteins', 0) > sequences_data[best_idx].get('Num_Source_Proteins', 0):
                    best_idx = idx
                    best_method = method
        
        representative = sequences_data[best_idx]
        
        # Aggregate ALL Source_Proteins
        all_proteins = []
        for epitope in cluster_epitopes:
            proteins_str = str(epitope.get('Source_Proteins', ''))
            if proteins_str and proteins_str != 'nan':
                proteins = [p.strip() for p in proteins_str.split(';')]
                all_proteins.extend(proteins)
        
        all_proteins = list(set(all_proteins))
        
        # Parse species
        species_set = set()
        for protein in all_proteins:
            species = parse_species_from_protein_id(protein)
            if species:
                species_set.add(species)
        
        species_list = sorted(list(species_set))
        num_species = len(species_list)
        
        # Aggregate metadata
        all_methods = set()
        all_epitope_ids = []
        all_data_sources = set()
        
        for epitope in cluster_epitopes:
            method = str(epitope.get('Method', ''))
            if method and method != 'nan':
                all_methods.add(method)
            
            epitope_id = epitope.get('Epitope_ID', '')
            if epitope_id:
                all_epitope_ids.append(epitope_id)
            
            data_source = str(epitope.get('Data_Source', ''))
            if data_source and data_source != 'nan':
                all_data_sources.add(data_source)
        
        # Save cluster summary
        cluster_results.append({
            'Cluster_ID': f"Cluster_{cluster_id}",
            'Cluster_Size': len(cluster_indices),
            'Representative_Sequence': representative['Epitope_Sequence'],
            'Representative_ID': representative.get('Epitope_ID', ''),
            'Representative_Length': len(representative['Epitope_Sequence']),
            'Representative_Method': representative.get('Method', ''),
            'Methods_All': '; '.join(sorted(all_methods)),
            'Has_IEDB': any('IEDB' in m for m in all_methods),
            'Num_Species': num_species,
            'Species_List': '; '.join(species_list),
            'Num_Source_Proteins': len(all_proteins),
            'Data_Sources': '; '.join(sorted(all_data_sources)),
            'Max_Human_Identity': representative.get('Max_Human_Identity', 0),
            'Max_Human_Coverage': representative.get('Max_Human_Coverage', 0),
            'All_Epitope_IDs': '; '.join(all_epitope_ids),
            'All_Sequences': '; '.join([e['Epitope_Sequence'] for e in cluster_epitopes])
        })
        
        # Detailed view
        for idx in cluster_indices:
            epitope = sequences_data[idx]
            detailed_clusters.append({
                'Cluster_ID': f"Cluster_{cluster_id}",
                'Is_Representative': (idx == best_idx),
                'Epitope_ID': epitope.get('Epitope_ID', ''),
                'Sequence': epitope['Epitope_Sequence'],
                'Length': len(epitope['Epitope_Sequence']),
                'Method': epitope.get('Method', ''),
                'Num_Source_Proteins': epitope.get('Num_Source_Proteins', 0),
                'Data_Source': epitope.get('Data_Source', ''),
                'Max_Human_Identity': epitope.get('Max_Human_Identity', 0)
            })
    
    print()
    
    # Create DataFrames
    cluster_summary_df = pd.DataFrame(cluster_results)
    detailed_df = pd.DataFrame(detailed_clusters)
    
    # Sort by species breadth
    cluster_summary_df = cluster_summary_df.sort_values('Num_Species', ascending=False)
    
    # Statistics
    print("\n" + "="*80)
    print("CLUSTER STATISTICS")
    print("="*80)
    
    print(f"\nTotal clusters: {len(cluster_summary_df):,}")
    print(f"Singleton clusters: {(cluster_summary_df['Cluster_Size'] == 1).sum():,}")
    print(f"Multi-member clusters: {(cluster_summary_df['Cluster_Size'] > 1).sum():,}")
    
    print(f"\nSpecies distribution:")
    print(f"  Mean: {cluster_summary_df['Num_Species'].mean():.1f}")
    print(f"  Median: {cluster_summary_df['Num_Species'].median():.0f}")
    print(f"  Max: {cluster_summary_df['Num_Species'].max()}")
    
    print(f"\nMethod distribution (representatives):")
    method_counts = cluster_summary_df['Representative_Method'].value_counts()
    for method, count in method_counts.head(10).items():
        pct = (count / len(cluster_summary_df)) * 100
        print(f"  {method:30} {count:6,} ({pct:5.1f}%)")
    
    # Top clusters
    print("\n" + "="*80)
    print("TOP 10 CLUSTERS BY SPECIES BREADTH")
    print("="*80)
    
    for _, row in cluster_summary_df.head(10).iterrows():
        print(f"\n{row['Cluster_ID']} ({row['Cluster_Size']} members)")
        print(f"  {row['Representative_Sequence']}")
        print(f"  Species: {row['Num_Species']} - {row['Species_List'][:100]}...")
        print(f"  Method: {row['Representative_Method']}")
    
    # Save
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    summary_file = f"{output_dir}/epitope_clusters_summary.txt"
    detailed_file = f"{output_dir}/epitope_clusters_detailed.txt"
    
    cluster_summary_df.to_csv(summary_file, sep='\t', index=False)
    detailed_df.to_csv(detailed_file, sep='\t', index=False)
    
    print(f"\n✓ Saved: {summary_file}")
    print(f"✓ Saved: {detailed_file}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    os.makedirs(output_dir, exist_ok=True)
    
    cluster_epitopes_fast(input_file, output_dir, log_dir, similarity_threshold=80)