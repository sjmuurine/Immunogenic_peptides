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
    """Parse species from protein ID - handles ANATOMIA and LITERATURE formats."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    # Species patterns (all 17 species)
    species_patterns = {
        'fragilis': 'B_fragilis', 'nucleatum': 'F_nucleatum', 'muciniphila': 'A_muciniphila',
        'prausnitzii': 'F_prausnitzii', 'reuteri': 'L_reuteri', 'cloacae': 'E_cloacae',
        'faecalis': 'E_faecalis', 'salivarius': 'S_salivarius', 'longum': 'B_longum',
        'plantarum': 'L_plantarum', 'casei': 'L_casei', 'sanguinis': 'T_sanguinis',
        'magna': 'V_magna', 'russellii': 'P_russellii', 'freudenreichii': 'P_freudenreichii',
        'fermentum': 'L_fermentum', 'adiacens': 'G_adiacens', 'gasseri': 'L_gasseri'
    }
    
    # LITERATURE format
    if 'LITERATURE' in protein_id:
        for pattern, species in species_patterns.items():
            if pattern.lower() in protein_id.lower():
                return species
        return None
    
    # ANATOMIA format: "ProteinID Species_OG_Conservation_Localization"
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

def cluster_sequences(sequences_data, similarity_threshold=80):
    """Cluster sequences by similarity."""
    n = len(sequences_data)
    clusters = []
    assigned = set()
    
    print(f"\n  Clustering {n} sequences with {similarity_threshold}% threshold...")
    
    for i in range(n):
        if i in assigned:
            continue
        
        if (i + 1) % 1000 == 0:
            print(f"    Processed {i+1}/{n} sequences...", end='\r')
        
        # Start new cluster
        cluster = [i]
        assigned.add(i)
        seq_i = sequences_data[i]['Epitope_Sequence']
        
        # Find similar sequences
        for j in range(i + 1, n):
            if j in assigned:
                continue
            
            seq_j = sequences_data[j]['Epitope_Sequence']
            
            # Calculate similarity
            similarity = calculate_simple_overlap_similarity(seq_i, seq_j)
            
            if similarity >= similarity_threshold:
                cluster.append(j)
                assigned.add(j)
        
        clusters.append(cluster)
    
    print()
    return clusters

def cluster_epitopes(input_file, output_dir, log_dir, similarity_threshold=80):
    """
    Cluster epitopes and aggregate metadata with corrected species parsing.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_06_cluster_epitopes")
    sys.stdout = logger
    
    print("="*80)
    print("EPITOPE CLUSTERING WITH SPECIES AGGREGATION")
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
    
    clusters = cluster_sequences(sequences_data, similarity_threshold)
    
    singleton_clusters = [c for c in clusters if len(c) == 1]
    multi_member_clusters = [c for c in clusters if len(c) > 1]
    
    print(f"\n  Total clusters: {len(clusters)}")
    print(f"  Singleton clusters: {len(singleton_clusters)}")
    print(f"  Multi-member clusters: {len(multi_member_clusters)}")
    
    if multi_member_clusters:
        sizes = [len(c) for c in multi_member_clusters]
        print(f"  Largest cluster: {max(sizes)} members")
        print(f"  Average multi-cluster size: {np.mean(sizes):.1f} members")
    
    # Analyze clusters and aggregate metadata
    print("\n" + "="*80)
    print("AGGREGATING METADATA AND PARSING SPECIES")
    print("="*80)
    
    cluster_results = []
    detailed_clusters = []
    
    for cluster_id, cluster_indices in enumerate(clusters, 1):
        if cluster_id % 1000 == 0:
            print(f"  Processing cluster {cluster_id}/{len(clusters)}...", end='\r')
        
        cluster_epitopes = [sequences_data[i] for i in cluster_indices]
        
        # Choose representative (IEDB priority, then by priority score)
        best_idx = cluster_indices[0]
        best_method = sequences_data[best_idx].get('Method', '')
        best_has_iedb = 'IEDB' in best_method
        best_priority = 0  # Will calculate later with species
        
        for idx in cluster_indices[1:]:
            method = sequences_data[idx].get('Method', '')
            has_iedb = 'IEDB' in method
            
            # Prioritize IEDB
            if has_iedb and not best_has_iedb:
                best_idx = idx
                best_method = method
                best_has_iedb = True
            elif has_iedb == best_has_iedb:
                # Both have or don't have IEDB, compare by num_source_proteins
                if sequences_data[idx].get('Num_Source_Proteins', 0) > sequences_data[best_idx].get('Num_Source_Proteins', 0):
                    best_idx = idx
                    best_method = method
        
        representative = sequences_data[best_idx]
        
        # Aggregate ALL Source_Proteins from cluster
        all_proteins = []
        for epitope in cluster_epitopes:
            proteins_str = str(epitope.get('Source_Proteins', ''))
            if proteins_str and proteins_str != 'nan':
                proteins = [p.strip() for p in proteins_str.split(';')]
                all_proteins.extend(proteins)
        
        # Remove duplicates
        all_proteins = list(set(all_proteins))
        
        # Parse species from aggregated proteins
        species_set = set()
        for protein in all_proteins:
            species = parse_species_from_protein_id(protein)
            if species:
                species_set.add(species)
        
        species_list = sorted(list(species_set))
        num_species = len(species_list)
        
        # Aggregate other metadata
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
        
        # Cluster summary
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
                'Max_Human_Identity': epitope.get('Max_Human_Identity', 0),
                'Max_Human_Coverage': epitope.get('Max_Human_Coverage', 0)
            })
    
    print()
    
    # Create DataFrames
    cluster_summary_df = pd.DataFrame(cluster_results)
    detailed_df = pd.DataFrame(detailed_clusters)
    
    # Sort by species breadth (for now, will re-rank later with full priority score)
    cluster_summary_df = cluster_summary_df.sort_values('Num_Species', ascending=False)
    
    # Statistics
    print("\n" + "="*80)
    print("CLUSTER STATISTICS")
    print("="*80)
    
    print(f"\nTotal clusters: {len(cluster_summary_df)}")
    print(f"Singleton clusters: {(cluster_summary_df['Cluster_Size'] == 1).sum()}")
    print(f"Multi-member clusters: {(cluster_summary_df['Cluster_Size'] > 1).sum()}")
    
    print(f"\nSpecies distribution:")
    print(f"  Mean species per cluster: {cluster_summary_df['Num_Species'].mean():.1f}")
    print(f"  Median species per cluster: {cluster_summary_df['Num_Species'].median():.0f}")
    print(f"  Max species per cluster: {cluster_summary_df['Num_Species'].max()}")
    
    print(f"\nMethod distribution:")
    method_counts = cluster_summary_df['Representative_Method'].value_counts()
    for method, count in method_counts.head(10).items():
        pct = (count / len(cluster_summary_df)) * 100
        print(f"  {method:30} {count:6} ({pct:5.1f}%)")
    
    # Show top clusters by species
    print("\n" + "="*80)
    print("TOP 10 CLUSTERS BY SPECIES BREADTH")
    print("="*80)
    
    top_clusters = cluster_summary_df.head(10)
    for _, row in top_clusters.iterrows():
        print(f"\n{row['Cluster_ID']} ({row['Cluster_Size']} members)")
        print(f"  Representative: {row['Representative_Sequence']}")
        print(f"  Species: {row['Num_Species']} - {row['Species_List']}")
        print(f"  Method: {row['Representative_Method']}")
        print(f"  IEDB: {row['Has_IEDB']}")
    
    # Save results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    summary_file = f"{output_dir}/epitope_clusters_summary.txt"
    detailed_file = f"{output_dir}/epitope_clusters_detailed.txt"
    
    cluster_summary_df.to_csv(summary_file, sep='\t', index=False)
    detailed_df.to_csv(detailed_file, sep='\t', index=False)
    
    print(f"\n✓ Saved summary: {summary_file}")
    print(f"✓ Saved detailed: {detailed_file}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("CLUSTERING COMPLETE")
    print("="*80)
    print("\nNext steps:")
    print("  1. Calculate priority scores for clusters")
    print("  2. Create tiers based on priority")
    print("  3. Select TOP candidates")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    os.makedirs(output_dir, exist_ok=True)
    
    cluster_epitopes(input_file, output_dir, log_dir, similarity_threshold=80)