import pandas as pd
import sys
import os
from datetime import datetime
import re

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

def extract_protein_id(protein_string):
    """
    Extract WP_XXXXXX.X from 'WP_XXXXXX.X Species_OG_Conservation_Localization'
    Returns None for LITERATURE entries.
    """
    protein_string = protein_string.strip()
    
    if protein_string.startswith('LITERATURE'):
        return None  # Skip LITERATURE entries
    elif protein_string.startswith('WP_'):
        # Extract just the WP_XXXXXX.X part (before space)
        return protein_string.split()[0]
    else:
        # For plain protein IDs, return as-is
        return protein_string

def merge_intensity_data_fixed(cluster_dir, intensity_analysis_dir, log_dir, anatomia_data_file):
    """
    Merge ANATOMIA intensity data with prioritized epitope clusters.
    Fixed version with proper protein ID extraction.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_24_merge_intensity_fixed")
    sys.stdout = logger
    
    print("="*80)
    print("MERGING ANATOMIA INTENSITY DATA (FIXED VERSION)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(intensity_analysis_dir, exist_ok=True)
    
    # Load ANATOMIA intensity data
    print(f"\nLoading ANATOMIA intensity data: {anatomia_data_file}")
    
    if os.path.exists(anatomia_data_file):
        anatomia_df = pd.read_csv(anatomia_data_file, sep='\t')
        print(f"✓ Loaded ANATOMIA data: {len(anatomia_df)} rows")
        print(f"  Columns: {list(anatomia_df.columns)}")
        print(f"  Sample protein IDs: {list(anatomia_df['Protein_ID'].head(3))}")
    else:
        print(f"❌ ANATOMIA data file not found: {anatomia_data_file}")
        logger.close()
        return
    
    # Load prioritized cluster data
    cluster_files = {
        'ALL': f"{cluster_dir}/epitope_clusters_prioritized_multi_species.txt",
        'TIER1': f"{cluster_dir}/TIER1_clusters_multi_species.txt",
        'TIER2': f"{cluster_dir}/TIER2_clusters_multi_species.txt"
    }
    
    print(f"\n" + "="*80)
    print("PROCESSING CLUSTER FILES")
    print("="*80)
    
    for file_type, file_path in cluster_files.items():
        if os.path.exists(file_path):
            print(f"\nProcessing {file_type}...")
            
            # Load cluster data
            clusters_df = pd.read_csv(file_path, sep='\t')
            print(f"  Loaded: {len(clusters_df)} clusters")
            
            # Create intensity summary per cluster
            cluster_intensities = []
            
            for cluster_idx, row in clusters_df.iterrows():
                cluster_id = row['Cluster_ID']
                proteins_str = str(row.get('Source_Proteins_All', ''))
                
                if proteins_str and proteins_str != 'nan':
                    # Split by semicolon and extract protein IDs
                    proteins_full = [p.strip() for p in proteins_str.split(';')]
                    protein_ids = []
                    literature_count = 0
                    
                    for protein_full in proteins_full:
                        protein_id = extract_protein_id(protein_full)
                        if protein_id:
                            protein_ids.append(protein_id)
                        elif protein_full.startswith('LITERATURE'):
                            literature_count += 1
                    
                    print(f"    {cluster_id}: {len(proteins_full)} total, {len(protein_ids)} ANATOMIA, {literature_count} LITERATURE")
                    
                    # Find matching ANATOMIA entries
                    matched_entries = []
                    matched_protein_ids = []
                    
                    for protein_id in protein_ids:
                        protein_matches = anatomia_df[anatomia_df['Protein_ID'] == protein_id]
                        if len(protein_matches) > 0:
                            matched_entries.extend(protein_matches.to_dict('records'))
                            matched_protein_ids.append(protein_id)
                    
                    # Calculate intensity statistics
                    if matched_entries:
                        # Find intensity columns
                        intensity_cols = [col for col in anatomia_df.columns 
                                        if 'intensity' in col.lower() or 'abundance' in col.lower()]
                        
                        intensity_stats = {
                            'Num_Total_Proteins': len(proteins_full),
                            'Num_ANATOMIA_Proteins': len(protein_ids),
                            'Num_LITERATURE_Proteins': literature_count,
                            'Num_Matched_ANATOMIA': len(matched_entries),
                            'Has_ANATOMIA_Data': True,
                            'Matched_Protein_IDs': '; '.join(matched_protein_ids[:10])  # Limit for readability
                        }
                        
                        # Add intensity statistics for each intensity column
                        for col in intensity_cols:
                            values = []
                            for entry in matched_entries:
                                val = entry.get(col)
                                if pd.notna(val) and val != '' and val != 0:
                                    try:
                                        values.append(float(val))
                                    except (ValueError, TypeError):
                                        pass
                            
                            if values:
                                intensity_stats[f'{col}_Mean'] = sum(values) / len(values)
                                intensity_stats[f'{col}_Max'] = max(values)
                                intensity_stats[f'{col}_Min'] = min(values)
                                intensity_stats[f'{col}_NonZero_Count'] = len(values)
                            else:
                                intensity_stats[f'{col}_Mean'] = 0
                                intensity_stats[f'{col}_Max'] = 0
                                intensity_stats[f'{col}_Min'] = 0
                                intensity_stats[f'{col}_NonZero_Count'] = 0
                        
                        print(f"      → Matched {len(matched_entries)} ANATOMIA proteins")
                    else:
                        intensity_stats = {
                            'Num_Total_Proteins': len(proteins_full),
                            'Num_ANATOMIA_Proteins': len(protein_ids),
                            'Num_LITERATURE_Proteins': literature_count,
                            'Num_Matched_ANATOMIA': 0,
                            'Has_ANATOMIA_Data': False,
                            'Matched_Protein_IDs': ''
                        }
                        print(f"      → No ANATOMIA matches found")
                
                else:
                    intensity_stats = {
                        'Num_Total_Proteins': 0,
                        'Num_ANATOMIA_Proteins': 0,
                        'Num_LITERATURE_Proteins': 0,
                        'Num_Matched_ANATOMIA': 0,
                        'Has_ANATOMIA_Data': False,
                        'Matched_Protein_IDs': ''
                    }
                
                cluster_intensities.append(intensity_stats)
            
            # Add intensity data to clusters
            intensity_df = pd.DataFrame(cluster_intensities)
            clusters_with_intensity = pd.concat([clusters_df, intensity_df], axis=1)
            
            # Save merged data
            output_file = f"{intensity_analysis_dir}/{file_type}_clusters_with_intensity_fixed.txt"
            clusters_with_intensity.to_csv(output_file, sep='\t', index=False)
            
            print(f"  ✓ Saved: {output_file}")
            
            # Statistics
            with_data = (clusters_with_intensity['Has_ANATOMIA_Data'] == True).sum()
            total_anatomia = clusters_with_intensity['Num_ANATOMIA_Proteins'].sum()
            matched_anatomia = clusters_with_intensity['Num_Matched_ANATOMIA'].sum()
            
            print(f"  Clusters with ANATOMIA data: {with_data}/{len(clusters_df)} ({with_data/len(clusters_df)*100:.1f}%)")
            if total_anatomia > 0:
                print(f"  ANATOMIA protein match rate: {matched_anatomia}/{total_anatomia} ({matched_anatomia/total_anatomia*100:.1f}%)")
            
        else:
            print(f"\n❌ File not found: {file_path}")
    
    print(f"\n" + "="*80)
    print("INTENSITY ANALYSIS SUMMARY")
    print("="*80)
    
    # Load TIER1 with intensity for detailed summary
    tier1_intensity_file = f"{intensity_analysis_dir}/TIER1_clusters_with_intensity_fixed.txt"
    if os.path.exists(tier1_intensity_file):
        tier1_df = pd.read_csv(tier1_intensity_file, sep='\t')
        
        print(f"\nTIER1 detailed analysis:")
        print(f"  Total TIER1 clusters: {len(tier1_df)}")
        
        with_data = (tier1_df['Has_ANATOMIA_Data'] == True).sum()
        print(f"  With ANATOMIA matches: {with_data} ({with_data/len(tier1_df)*100:.1f}%)")
        
        if with_data > 0:
            # Show intensity column statistics
            intensity_cols = [col for col in tier1_df.columns if col.endswith('_Mean') and 'Intensity' in col]
            
            print(f"\n  Intensity statistics for TIER1 with data:")
            for col in intensity_cols:
                values = tier1_df[tier1_df['Has_ANATOMIA_Data'] == True][col]
                non_zero = values[values > 0]
                if len(non_zero) > 0:
                    print(f"    {col}: {len(non_zero)} non-zero, mean={non_zero.mean():.6f}, max={non_zero.max():.6f}")
            
            # Show top 10 by intensity
            if len(intensity_cols) > 0:
                main_intensity_col = intensity_cols[0]  # Use first intensity column
                top_intensity = tier1_df[tier1_df['Has_ANATOMIA_Data'] == True].nlargest(10, main_intensity_col)
                
                print(f"\n  Top 10 TIER1 by {main_intensity_col}:")
                print(f"  {'Sequence':<35} {'Sp':<3} {'ID%':<6} {'Intensity':<10}")
                print("  " + "-" * 70)
                for _, row in top_intensity.iterrows():
                    seq = row['Representative_Sequence'][:33]
                    species = row['Num_Species']
                    identity = row['Max_Human_Identity']
                    intensity = row[main_intensity_col]
                    print(f"  {seq:<35} {species:<3} {identity:<6.1f} {intensity:<10.6f}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("FIXED INTENSITY DATA MERGING COMPLETE")
    print("="*80)
    print("\nOutput files:")
    print(f"  {intensity_analysis_dir}/ALL_clusters_with_intensity_fixed.txt")
    print(f"  {intensity_analysis_dir}/TIER1_clusters_with_intensity_fixed.txt")  
    print(f"  {intensity_analysis_dir}/TIER2_clusters_with_intensity_fixed.txt")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    # Paths
    cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    anatomia_data_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'
    
    merge_intensity_data_fixed(cluster_dir, intensity_analysis_dir, log_dir, anatomia_data_file)