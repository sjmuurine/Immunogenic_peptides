#!/bin/bash
#SBATCH --job-name=cluster_epitopes
#SBATCH --account=project_2009813
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_epitopes_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_epitopes_%j.err

echo "=================================================="
echo "Fast Epitope Clustering"
echo "=================================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load required modules
echo "Loading modules..."
module load gcc/11.3.0
module load biopythontools
echo "Modules loaded"
echo ""

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os
from datetime import datetime

def parse_species_from_protein_id(protein_id):
    """Parse species from protein ID."""
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

def cluster_sequences_memory_efficient(sequences_data, similarity_threshold=80):
    """Memory-efficient clustering - stream comparisons instead of storing."""
    
    n = len(sequences_data)
    print(f"\n  Clustering {n} sequences with {similarity_threshold}% threshold...")
    
    # Group by length
    print("  Step 1: Grouping by length...")
    length_groups = defaultdict(list)
    for i, data in enumerate(sequences_data):
        length = len(data['Epitope_Sequence'])
        length_groups[length].append(i)
    
    print(f"    Found {len(length_groups)} unique lengths")
    
    # Cluster - stream comparisons
    print("  Step 2: Clustering (memory-efficient streaming)...")
    
    clusters = []
    assigned = set()
    comparisons = 0
    
    for i in range(n):
        if i in assigned:
            continue
        
        if len(clusters) % 5000 == 0 and len(clusters) > 0:
            print(f"    Progress: {len(clusters):,} clusters, {comparisons:,} comparisons", end='\r')
        
        # Start new cluster
        cluster = [i]
        assigned.add(i)
        seq_i = sequences_data[i]['Epitope_Sequence']
        len_i = len(seq_i)
        
        # Only compare with nearby lengths (stream, don't store)
        for nearby_len in range(len_i - 2, len_i + 3):
            if nearby_len not in length_groups:
                continue
            
            for j in length_groups[nearby_len]:
                if j <= i or j in assigned:
                    continue
                
                seq_j = sequences_data[j]['Epitope_Sequence']
                
                # Quick substring check
                if seq_i in seq_j or seq_j in seq_i:
                    cluster.append(j)
                    assigned.add(j)
                    comparisons += 1
                    continue
                
                # Length pre-filter
                len_diff = abs(len_i - len(seq_j))
                if len_diff > max(len_i, len(seq_j)) * 0.2:
                    comparisons += 1
                    continue
                
                # Full similarity calc
                similarity = calculate_simple_overlap_similarity(seq_i, seq_j)
                comparisons += 1
                
                if similarity >= similarity_threshold:
                    cluster.append(j)
                    assigned.add(j)
        
        clusters.append(cluster)
    
    print()
    print(f"  Total comparisons: {comparisons:,}")
    print(f"  Clusters created: {len(clusters):,}")
    
    return clusters

# Main execution
print("="*80)
print("FAST EPITOPE CLUSTERING")
print("="*80)
print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Load data
input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
print(f"\nLoading: {input_file}")

epitopes_df = pd.read_csv(input_file, sep='\t')
print(f"Total epitopes: {len(epitopes_df):,}")

# Filter by length
epitopes_df = epitopes_df[
    (epitopes_df['Epitope_Sequence'].str.len() >= 11) & 
    (epitopes_df['Epitope_Sequence'].str.len() <= 24)
]
print(f"After length filter (11-24 aa): {len(epitopes_df):,}")

sequences_data = epitopes_df.to_dict('records')

# Cluster
print("\n" + "="*80)
print("CLUSTERING")
print("="*80)

clusters = cluster_sequences_memory_efficient(sequences_data, similarity_threshold=80)

print(f"\nClustering complete!")
print(f"  Total clusters: {len(clusters):,}")
print(f"  Singletons: {sum(1 for c in clusters if len(c) == 1):,}")
print(f"  Multi-member: {sum(1 for c in clusters if len(c) > 1):,}")

# Aggregate metadata
print("\n" + "="*80)
print("AGGREGATING METADATA")
print("="*80)

cluster_results = []

for cluster_id, cluster_indices in enumerate(clusters, 1):
    if cluster_id % 10000 == 0:
        print(f"  Processing {cluster_id:,}/{len(clusters):,}...", end='\r')
    
    cluster_epitopes = [sequences_data[i] for i in cluster_indices]
    
    # Choose representative (IEDB priority)
    best_idx = cluster_indices[0]
    for idx in cluster_indices[1:]:
        if 'IEDB' in sequences_data[idx].get('Method', '') and 'IEDB' not in sequences_data[best_idx].get('Method', ''):
            best_idx = idx
    
    representative = sequences_data[best_idx]
    
    # Aggregate proteins
    all_proteins = []
    for epitope in cluster_epitopes:
        proteins_str = str(epitope.get('Source_Proteins', ''))
        if proteins_str and proteins_str != 'nan':
            all_proteins.extend([p.strip() for p in proteins_str.split(';')])
    
    all_proteins = list(set(all_proteins))
    
    # Parse species
    species_set = set()
    for protein in all_proteins:
        species = parse_species_from_protein_id(protein)
        if species:
            species_set.add(species)
    
    # Aggregate methods
    all_methods = set()
    for epitope in cluster_epitopes:
        method = str(epitope.get('Method', ''))
        if method and method != 'nan':
            all_methods.add(method)
    
    cluster_results.append({
        'Cluster_ID': f"Cluster_{cluster_id}",
        'Cluster_Size': len(cluster_indices),
        'Representative_Sequence': representative['Epitope_Sequence'],
        'Representative_ID': representative.get('Epitope_ID', ''),
        'Representative_Method': representative.get('Method', ''),
        'Methods_All': '; '.join(sorted(all_methods)),
        'Has_IEDB': any('IEDB' in m for m in all_methods),
        'Num_Species': len(species_set),
        'Species_List': '; '.join(sorted(species_set)),
        'Num_Source_Proteins': len(all_proteins),
        'Max_Human_Identity': representative.get('Max_Human_Identity', 0),
        'Max_Human_Coverage': representative.get('Max_Human_Coverage', 0)
    })

print()

# Save
cluster_df = pd.DataFrame(cluster_results)
cluster_df = cluster_df.sort_values('Num_Species', ascending=False)

output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
os.makedirs(output_dir, exist_ok=True)

output_file = f"{output_dir}/epitope_clusters_summary.txt"
cluster_df.to_csv(output_file, sep='\t', index=False)

print(f"\n✓ Saved: {output_file}")
print(f"\nTop 5 clusters by species:")
for _, row in cluster_df.head(5).iterrows():
    print(f"  {row['Representative_Sequence'][:30]:30} - {row['Num_Species']} species")

print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

PYTHON_SCRIPT

echo ""
echo "=================================================="
echo "Clustering Complete"
echo "End time: $(date)"
echo "=================================================="