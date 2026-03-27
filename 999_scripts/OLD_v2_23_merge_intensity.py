import pandas as pd
import sys
import os
from datetime import datetime
import shutil

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

def merge_intensity_data(cluster_dir, intensity_analysis_dir, log_dir, anatomia_data_file):
    """
    Merge ANATOMIA intensity data with prioritized epitope clusters.
    """
    
    # Set up logging
    logger = Logger(log_dir, "v2_23_merge_intensity")
    sys.stdout = logger
    
    print("="*80)
    print("MERGING ANATOMIA INTENSITY DATA WITH EPITOPE CLUSTERS")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(intensity_analysis_dir, exist_ok=True)
    
    # Load ANATOMIA intensity data
    print(f"\nLoading ANATOMIA intensity data...")
    print(f"Expected file: {anatomia_data_file}")
    
    if os.path.exists(anatomia_data_file):
        # Try different formats (Excel, CSV, TSV)
        if anatomia_data_file.endswith('.xlsx') or anatomia_data_file.endswith('.xls'):
            anatomia_df = pd.read_excel(anatomia_data_file)
        elif anatomia_data_file.endswith('.csv'):
            anatomia_df = pd.read_csv(anatomia_data_file)
        else:
            anatomia_df = pd.read_csv(anatomia_data_file, sep='\t')
        
        print(f"✓ Loaded ANATOMIA data: {len(anatomia_df)} rows")
        print(f"  Columns: {list(anatomia_df.columns)}")
    else:
        print(f"❌ ANATOMIA data file not found!")
        print(f"Please provide the path to your ANATOMIA intensity data file")
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
            
            # Extract all source proteins from clusters
            all_proteins = []
            for idx, row in clusters_df.iterrows():
                proteins_str = str(row.get('Source_Proteins_All', ''))
                if proteins_str and proteins_str != 'nan':
                    proteins = [p.strip() for p in proteins_str.split(';')]
                    all_proteins.extend([(idx, p) for p in proteins])
            
            print(f"  Extracted {len(all_proteins)} protein-cluster pairs")
            
            # Create protein mapping
            protein_to_cluster = {}
            for cluster_idx, protein_id in all_proteins:
                if protein_id not in protein_to_cluster:
                    protein_to_cluster[protein_id] = []
                protein_to_cluster[protein_id].append(cluster_idx)
            
            # Merge with ANATOMIA data
            print(f"  Merging with ANATOMIA intensity data...")
            
            # Create intensity summary per cluster
            cluster_intensities = []
            
            for cluster_idx, row in clusters_df.iterrows():
                cluster_id = row['Cluster_ID']
                proteins_str = str(row.get('Source_Proteins_All', ''))
                
                if proteins_str and proteins_str != 'nan':
                    proteins = [p.strip() for p in proteins_str.split(';')]
                    
                    # Find matching ANATOMIA entries
                    matched_intensities = []
                    for protein_id in proteins:
                        # Try exact match first
                        protein_matches = anatomia_df[anatomia_df['Protein_ID'] == protein_id]
                        if len(protein_matches) == 0:
                            # Try partial match (in case of format differences)
                            protein_matches = anatomia_df[anatomia_df['Protein_ID'].str.contains(protein_id, na=False)]
                        
                        if len(protein_matches) > 0:
                            matched_intensities.extend(protein_matches.to_dict('records'))
                    
                    # Calculate intensity statistics
                    if matched_intensities:
                        intensity_cols = [col for col in anatomia_df.columns if 'intensity' in col.lower() or 'abundance' in col.lower()]
                        
                        intensity_stats = {
                            'Num_ANATOMIA_Proteins': len(matched_intensities),
                            'Has_ANATOMIA_Data': True
                        }
                        
                        # Add intensity statistics if intensity columns exist
                        for col in intensity_cols[:5]:  # Limit to first 5 intensity columns
                            values = [m.get(col, 0) for m in matched_intensities if pd.notna(m.get(col, None))]
                            if values:
                                intensity_stats[f'{col}_Mean'] = sum(values) / len(values)
                                intensity_stats[f'{col}_Max'] = max(values)
                                intensity_stats[f'{col}_Min'] = min(values)
                    else:
                        intensity_stats = {
                            'Num_ANATOMIA_Proteins': 0,
                            'Has_ANATOMIA_Data': False
                        }
                
                else:
                    intensity_stats = {
                        'Num_ANATOMIA_Proteins': 0,
                        'Has_ANATOMIA_Data': False
                    }
                
                cluster_intensities.append(intensity_stats)
            
            # Add intensity data to clusters
            intensity_df = pd.DataFrame(cluster_intensities)
            clusters_with_intensity = pd.concat([clusters_df, intensity_df], axis=1)
            
            # Save merged data
            output_file = f"{intensity_analysis_dir}/{file_type}_clusters_with_intensity.txt"
            clusters_with_intensity.to_csv(output_file, sep='\t', index=False)
            
            print(f"  ✓ Saved: {output_file}")
            
            # Statistics
            with_data = (clusters_with_intensity['Has_ANATOMIA_Data'] == True).sum()
            print(f"  Clusters with ANATOMIA data: {with_data}/{len(clusters_df)} ({with_data/len(clusters_df)*100:.1f}%)")
            
        else:
            print(f"\n❌ File not found: {file_path}")
    
    print(f"\n" + "="*80)
    print("INTENSITY ANALYSIS SUMMARY")
    print("="*80)
    
    # Load TIER1 with intensity for summary
    tier1_intensity_file = f"{intensity_analysis_dir}/TIER1_clusters_with_intensity.txt"
    if os.path.exists(tier1_intensity_file):
        tier1_df = pd.read_csv(tier1_intensity_file, sep='\t')
        
        print(f"\nTIER1 intensity analysis:")
        print(f"  Total TIER1 clusters: {len(tier1_df)}")
        
        with_data = (tier1_df['Has_ANATOMIA_Data'] == True).sum()
        print(f"  With ANATOMIA data: {with_data} ({with_data/len(tier1_df)*100:.1f}%)")
        
        if with_data > 0:
            print(f"  Mean proteins per cluster: {tier1_df[tier1_df['Has_ANATOMIA_Data']]['Num_ANATOMIA_Proteins'].mean():.1f}")
            
            # Show top clusters by ANATOMIA data
            top_anatomia = tier1_df[tier1_df['Has_ANATOMIA_Data'] == True].head(10)
            print(f"\n  Top TIER1 clusters with ANATOMIA data:")
            print(f"  {'Sequence':<35} {'Sp':<3} {'ID%':<6} {'ANATOMIA'}")
            print("  " + "-" * 65)
            for _, row in top_anatomia.iterrows():
                seq = row['Representative_Sequence'][:33]
                species = row['Num_Species']
                identity = row['Max_Human_Identity']
                anatomia_count = row['Num_ANATOMIA_Proteins']
                print(f"  {seq:<35} {species:<3} {identity:<6.1f} {anatomia_count}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    print("INTENSITY DATA MERGING COMPLETE")
    print("="*80)
    print("\nOutput files:")
    print(f"  {intensity_analysis_dir}/ALL_clusters_with_intensity.txt")
    print(f"  {intensity_analysis_dir}/TIER1_clusters_with_intensity.txt")  
    print(f"  {intensity_analysis_dir}/TIER2_clusters_with_intensity.txt")
    print("\nNext steps:")
    print("  1. Review intensity data for TIER1 candidates")
    print("  2. Visualize intensity patterns")
    print("  3. Proceed with metaproteome validation")
    print("="*80)
    
    logger.close()

if __name__ == "__main__":
    # Paths
    cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    # You'll need to provide the path to your ANATOMIA intensity data
    anatomia_data_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'  # Update this path!
    
    print("="*80)
    print("ANATOMIA INTENSITY DATA MERGER")
    print("="*80)
    print("\n⚠️  IMPORTANT: Please update the ANATOMIA data file path in the script!")
    print(f"Current path: {anatomia_data_file}")
    print("Make sure this points to your ANATOMIA intensity/abundance data file.")
    print("\nPress Ctrl+C to cancel and update the path, or continue if path is correct.")
    
    try:
        merge_intensity_data(cluster_dir, intensity_analysis_dir, log_dir, anatomia_data_file)
    except KeyboardInterrupt:
        print("\n\nScript cancelled. Please update the ANATOMIA data file path and run again.")