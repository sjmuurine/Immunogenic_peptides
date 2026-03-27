import pandas as pd
import numpy as np
from collections import defaultdict

def calculate_simple_overlap_similarity(seq1, seq2):
    """
    Calculate similarity by trying all possible alignments.
    Handles frame shifts and overlaps.
    """
    if seq1 == seq2:
        return 100.0
    
    max_similarity = 0
    
    # Try all offsets
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
            # Require at least 80% of shorter sequence to overlap
            min_len = min(len(seq1), len(seq2))
            if overlap_length >= min_len * 0.8:
                similarity = (matches / overlap_length) * 100
                max_similarity = max(max_similarity, similarity)
    
    return max_similarity

def cluster_sequences(sequences_data, similarity_threshold=70):
    """
    Cluster sequences by similarity.
    """
    n = len(sequences_data)
    clusters = []
    assigned = set()
    
    print(f"  Clustering {n} sequences...")
    
    for i in range(n):
        if i in assigned:
            continue
        
        if (i + 1) % 50 == 0:
            print(f"    Processed {i+1}/{n} sequences...")
        
        # Start new cluster
        cluster = [i]
        assigned.add(i)
        seq_i = sequences_data[i]['Sequence']
        
        # Find similar sequences
        for j in range(i + 1, n):
            if j in assigned:
                continue
            
            seq_j = sequences_data[j]['Sequence']
            
            # Calculate similarity
            similarity = calculate_simple_overlap_similarity(seq_i, seq_j)
            
            if similarity >= similarity_threshold:
                cluster.append(j)
                assigned.add(j)
        
        clusters.append(cluster)
    
    return clusters

def analyze_epitope_clusters(input_dir, output_dir, similarity_threshold=75):
    """
    Cluster epitopes with improved similarity detection.
    """
    
    print("="*70)
    print("EPITOPE CLUSTERING ANALYSIS")
    print(f"Similarity threshold: {similarity_threshold}%")
    print("="*70)
    
    # Files to process
    files_to_cluster = {
        'TOP50': f"{input_dir}/TOP50_epitopes.txt",
        'TIER1': f"{input_dir}/TIER1_top_priority_epitopes.txt",
        'TIER2': f"{input_dir}/TIER2_high_priority_epitopes.txt",
        'COMBINED_TOP_TIERS': None
    }
    
    # Load protein mappings
    protein_mappings = {}
    for tier in ['TOP50', 'TIER1', 'TIER2']:
        mapping_file = f"{input_dir}/{tier}_protein_mapping.txt"
        try:
            protein_mappings[tier] = pd.read_csv(mapping_file, sep='\t')
        except:
            protein_mappings[tier] = pd.DataFrame()
    
    # Process each file
    for dataset_name, file_path in files_to_cluster.items():
        print(f"\n{'='*70}")
        print(f"Processing: {dataset_name}")
        print("="*70)
        
        # Load data
        if dataset_name == 'COMBINED_TOP_TIERS':
            dfs = []
            for tier in ['TOP50', 'TIER1', 'TIER2']:
                try:
                    df = pd.read_csv(files_to_cluster[tier], sep='\t')
                    dfs.append(df)
                except:
                    pass
            
            if not dfs:
                continue
            
            epitopes_df = pd.concat(dfs, ignore_index=True)
            epitopes_df = epitopes_df.drop_duplicates(subset=['Epitope_ID'])
        else:
            try:
                epitopes_df = pd.read_csv(file_path, sep='\t')
            except:
                continue
        
        print(f"Total epitopes: {len(epitopes_df)}")
        
        if len(epitopes_df) == 0:
            continue
        
        # Prepare data
        sequences_data = epitopes_df.to_dict('records')
        
        # Perform clustering
        clusters = cluster_sequences(sequences_data, similarity_threshold)
        
        multi_member_clusters = [c for c in clusters if len(c) > 1]
        print(f"\n  Found {len(clusters)} total clusters")
        print(f"  Multi-variant clusters: {len(multi_member_clusters)}")
        
        # Analyze clusters
        cluster_results = []
        detailed_clusters = []
        
        for cluster_id, cluster_indices in enumerate(clusters, 1):
            cluster_epitopes = [sequences_data[i] for i in cluster_indices]
            sequences = [e['Sequence'] for e in cluster_epitopes]
            epitope_ids = [e['Epitope_ID'] for e in cluster_epitopes]
            
            # Pick best representative
            best_idx = cluster_indices[0]
            best_priority = sequences_data[best_idx].get('Priority_Score', 0)
            for idx in cluster_indices[1:]:
                score = sequences_data[idx].get('Priority_Score', 0)
                if score > best_priority:
                    best_idx = idx
                    best_priority = score
            
            representative = sequences_data[best_idx]
            
            # Aggregate metadata
            total_hits = sum(e.get('Total_Hits', 0) for e in cluster_epitopes)
            anatomia_hits = sum(e.get('ANATOMIA_Hits', 0) for e in cluster_epitopes)
            literature_hits = sum(e.get('LITERATURE_Hits', 0) for e in cluster_epitopes)
            
            all_methods = set()
            for e in cluster_epitopes:
                methods = str(e.get('Source_Methods', '')).split(',')
                all_methods.update(methods)
            
            # Get protein mappings
            cluster_proteins = []
            for epitope_id in epitope_ids:
                for tier, mapping_df in protein_mappings.items():
                    if len(mapping_df) > 0:
                        matches = mapping_df[mapping_df['Epitope_ID'] == epitope_id]
                        if len(matches) > 0:
                            cluster_proteins.append(matches)
            
            if cluster_proteins:
                combined_proteins = pd.concat(cluster_proteins, ignore_index=True)
                combined_proteins = combined_proteins.drop_duplicates(subset=['Protein_ID'])
                
                unique_species = combined_proteins[
                    (combined_proteins['Species'] != 'Unknown') & 
                    (combined_proteins['Species'] != 'LITERATURE')
                ]['Species'].nunique()
                
                unique_ogs = combined_proteins[
                    (combined_proteins['Orthogroup'] != 'Unknown') & 
                    (combined_proteins['Orthogroup'] != 'N/A')
                ]['Orthogroup'].nunique()
                
                total_proteins = len(combined_proteins)
            else:
                unique_species = unique_ogs = total_proteins = 0
            
            cluster_results.append({
                'Cluster_ID': f"Cluster_{cluster_id}",
                'Cluster_Size': len(cluster_indices),
                'Representative_Sequence': representative['Sequence'],
                'Representative_ID': representative['Epitope_ID'],
                'Representative_Length': representative.get('Length', len(representative['Sequence'])),
                'Total_Hits_AllVariants': total_hits,
                'ANATOMIA_Hits_AllVariants': anatomia_hits,
                'LITERATURE_Hits_AllVariants': literature_hits,
                'Total_Proteins_AllVariants': total_proteins,
                'Unique_Species_AllVariants': unique_species,
                'Unique_OGs_AllVariants': unique_ogs,
                'Methods_Used': ','.join(sorted(all_methods)),
                'Has_HIGH_CONFIDENCE': 'HIGH_CONFIDENCE' in all_methods,
                'Has_IEDB': any('IEDB' in m for m in all_methods),
                'Has_Kmer_EXACT': 'Kmer_EXACT_only' in all_methods,
                'Representative_Priority_Score': best_priority,
                'All_Sequences': '; '.join(sequences),
                'Human_Identity': representative.get('Max_Human_Identity', 0)
            })
            
            for idx in cluster_indices:
                epitope = sequences_data[idx]
                detailed_clusters.append({
                    'Cluster_ID': f"Cluster_{cluster_id}",
                    'Is_Representative': (idx == best_idx),
                    'Epitope_ID': epitope['Epitope_ID'],
                    'Sequence': epitope['Sequence'],
                    'Length': epitope.get('Length', len(epitope['Sequence'])),
                    'Source_Methods': epitope.get('Source_Methods', ''),
                    'ANATOMIA_Hits': epitope.get('ANATOMIA_Hits', 0),
                    'LITERATURE_Hits': epitope.get('LITERATURE_Hits', 0),
                    'Total_Hits': epitope.get('Total_Hits', 0),
                    'Priority_Score': epitope.get('Priority_Score', 0),
                    'Human_Identity': epitope.get('Max_Human_Identity', 0)
                })
        
        # Sort and save
        cluster_summary_df = pd.DataFrame(cluster_results)
        cluster_summary_df = cluster_summary_df.sort_values('Representative_Priority_Score', ascending=False)
        detailed_df = pd.DataFrame(detailed_clusters)
        
        summary_file = f"{output_dir}/{dataset_name}_cluster_summary_v2.txt"
        detailed_file = f"{output_dir}/{dataset_name}_cluster_detailed_v2.txt"
        
        cluster_summary_df.to_csv(summary_file, sep='\t', index=False)
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        
        print(f"\n✓ Saved: {summary_file}")
        print(f"✓ Saved: {detailed_file}")
        
        print(f"\nCluster Statistics:")
        print(f"  Singleton clusters: {len(cluster_summary_df[cluster_summary_df['Cluster_Size'] == 1])}")
        print(f"  Multi-variant clusters: {len(cluster_summary_df[cluster_summary_df['Cluster_Size'] > 1])}")
        if len(cluster_summary_df) > 0:
            print(f"  Largest cluster size: {cluster_summary_df['Cluster_Size'].max()}")
        
        # Show multi-member clusters
        multi_clusters = cluster_summary_df[cluster_summary_df['Cluster_Size'] > 1].head(10)
        if len(multi_clusters) > 0:
            print(f"\nTop multi-variant clusters:")
            for _, row in multi_clusters.iterrows():
                print(f"\n  {row['Cluster_ID']} (size {row['Cluster_Size']}):")
                sequences = row['All_Sequences'].split('; ')
                for seq in sequences[:5]:  # Show first 5
                    print(f"    {seq}")
                if len(sequences) > 5:
                    print(f"    ... and {len(sequences)-5} more")
    
    print("\n" + "="*70)
    print("CLUSTERING COMPLETE")
    print("="*70)

# Run
input_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
analyze_epitope_clusters(input_dir, output_dir, similarity_threshold=75)