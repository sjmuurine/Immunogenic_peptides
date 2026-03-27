import pandas as pd
import sys
import os
from datetime import datetime
import shutil

def create_refined_tier1_for_metaproteome(intensity_analysis_dir, metaproteome_dir, log_dir):
    """
    Create refined TIER1 combining original TIER1 + high-abundance TIER2/3 candidates.
    """
    
    # Set up logging
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"{log_dir}/{timestamp}_create_refined_tier1.log"
    
    print("="*80)
    print("CREATING REFINED TIER1 FOR METAPROTEOME VALIDATION")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Log: {log_file}")
    
    # Create metaproteome directory
    os.makedirs(metaproteome_dir, exist_ok=True)
    
    # Load abundance data with high-abundance candidates
    abundance_file = f"{intensity_analysis_dir}/ALL_clusters_with_intensity_fixed.txt"
    print(f"\nLoading abundance data: {abundance_file}")
    
    df = pd.read_csv(abundance_file, sep='\t')
    anatomia_df = df[df['Has_ANATOMIA_Data'] == True].copy()
    
    # Calculate total abundance (same logic as before)
    intensity_file = f"{intensity_analysis_dir}/../001_data/anatomia_protein_intensities.txt"
    intensity_raw_df = pd.read_csv(intensity_file, sep='\t')
    protein_intensity_map = dict(zip(intensity_raw_df['Protein_ID'], 
                                   intensity_raw_df['Protein_Mean_Relative_Intensity']))
    
    def calculate_total_abundance(proteins_str):
        if pd.isna(proteins_str) or proteins_str == '':
            return 0.0
        total_abundance = 0.0
        for protein_full in proteins_str.split(';'):
            protein_full = protein_full.strip()
            if protein_full.startswith('WP_'):
                protein_id = protein_full.split()[0]
                if protein_id in protein_intensity_map:
                    abundance = protein_intensity_map[protein_id]
                    if pd.notna(abundance) and abundance > 0:
                        total_abundance += abundance
        return total_abundance
    
    anatomia_df['Total_Relative_Abundance'] = anatomia_df['Source_Proteins_All'].apply(calculate_total_abundance)
    abundance_df = anatomia_df[anatomia_df['Total_Relative_Abundance'] > 0].copy()
    
    print(f"Clusters with abundance data: {len(abundance_df)}")
    
    # Step 1: Get current TIER1
    current_tier1 = abundance_df[abundance_df['Tier'] == 'TIER1'].copy()
    print(f"\nCurrent TIER1: {len(current_tier1)} clusters")
    
    # Step 2: Find high-abundance TIER2/3 candidates with ≤75% human identity
    tier1_abundances = current_tier1['Total_Relative_Abundance']
    if len(tier1_abundances) > 0:
        abundance_threshold = tier1_abundances.quantile(0.75)
    else:
        abundance_threshold = abundance_df['Total_Relative_Abundance'].quantile(0.9)
    
    high_abundance_candidates = abundance_df[
        (abundance_df['Tier'].isin(['TIER2', 'TIER3'])) & 
        (abundance_df['Total_Relative_Abundance'] > abundance_threshold) &
        (abundance_df['Max_Human_Identity'] <= 75.0)  # ≤75% human identity
    ].copy()
    
    print(f"\nHigh-abundance TIER2/3 candidates (≤75% human identity):")
    print(f"  Abundance threshold: {abundance_threshold:.4f}")
    print(f"  Found: {len(high_abundance_candidates)} candidates")
    
    if len(high_abundance_candidates) > 0:
        print(f"\n  High-abundance candidates to promote:")
        print(f"  {'Sequence':<35} {'Tier':<6} {'Sp':<3} {'ID%':<6} {'Abundance':<12}")
        print("  " + "-" * 80)
        for _, row in high_abundance_candidates.sort_values('Total_Relative_Abundance', ascending=False).iterrows():
            seq = row['Representative_Sequence'][:33]
            tier = row['Tier']
            species = row['Num_Species']
            identity = row['Max_Human_Identity']
            abundance = row['Total_Relative_Abundance']
            print(f"  {seq:<35} {tier:<6} {species:<3} {identity:<6.1f} {abundance:<12.4f}")
    
    # Step 3: Combine into refined TIER1
    refined_tier1 = pd.concat([current_tier1, high_abundance_candidates], ignore_index=True)
    refined_tier1 = refined_tier1.sort_values('Priority_Score', ascending=False)
    refined_tier1['Refined_Tier'] = 'REFINED_TIER1'
    
    print(f"\n" + "="*80)
    print("REFINED TIER1 SUMMARY")
    print("="*80)
    print(f"Original TIER1: {len(current_tier1)}")
    print(f"Added from TIER2/3: {len(high_abundance_candidates)}")
    print(f"Refined TIER1 total: {len(refined_tier1)}")
    
    # Step 4: Save refined TIER1
    refined_tier1_file = f"{intensity_analysis_dir}/REFINED_TIER1_for_metaproteome.txt"
    refined_tier1.to_csv(refined_tier1_file, sep='\t', index=False)
    print(f"\n✓ Saved refined TIER1: {refined_tier1_file}")
    
    # Copy to metaproteome directory
    metaproteome_tier1_file = f"{metaproteome_dir}/REFINED_TIER1_clusters.txt"
    shutil.copy2(refined_tier1_file, metaproteome_tier1_file)
    print(f"✓ Copied to metaproteome dir: {metaproteome_tier1_file}")
    
    # Step 5: Create FASTA file for BLAST analysis
    fasta_file = f"{metaproteome_dir}/refined_tier1_epitopes.fasta"
    print(f"\nCreating FASTA file for BLAST: {fasta_file}")
    
    with open(fasta_file, 'w') as f:
        for _, row in refined_tier1.iterrows():
            cluster_id = row['Cluster_ID']
            sequence = row['Representative_Sequence']
            f.write(f">{cluster_id}\n{sequence}\n")
    
    print(f"✓ Created FASTA: {fasta_file} ({len(refined_tier1)} sequences)")
    
    # Step 6: Summary statistics
    print(f"\n" + "="*80)
    print("REFINED TIER1 ANALYSIS")
    print("="*80)
    
    print(f"\nHuman identity distribution:")
    for threshold in [60, 65, 70, 75, 80]:
        count = (refined_tier1['Max_Human_Identity'] <= threshold).sum()
        pct = count / len(refined_tier1) * 100
        print(f"  ≤{threshold}%: {count}/{len(refined_tier1)} ({pct:.1f}%)")
    
    print(f"\nSpecies distribution:")
    print(f"  Mean: {refined_tier1['Num_Species'].mean():.1f}")
    print(f"  Range: {refined_tier1['Num_Species'].min()}-{refined_tier1['Num_Species'].max()}")
    
    print(f"\nAbundance distribution:")
    print(f"  Mean: {refined_tier1['Total_Relative_Abundance'].mean():.4f}")
    print(f"  Range: {refined_tier1['Total_Relative_Abundance'].min():.4f}-{refined_tier1['Total_Relative_Abundance'].max():.4f}")
    
    print(f"\nMethod distribution:")
    method_counts = refined_tier1['Representative_Method'].value_counts()
    for method, count in method_counts.items():
        pct = count / len(refined_tier1) * 100
        print(f"  {method}: {count} ({pct:.1f}%)")
    
    print(f"\n" + "="*80)
    print("NEXT STEPS: METAPROTEOME VALIDATION")
    print("="*80)
    print(f"\nFiles ready for metaproteome analysis:")
    print(f"  Epitope clusters: {metaproteome_tier1_file}")
    print(f"  FASTA for BLAST: {fasta_file}")
    print(f"  Metaproteome data: /scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/metaproteome_sequences.fa")
    print(f"\nSuggested approaches:")
    print(f"  1. BLAST search (allows mismatches)")
    print(f"  2. Exact string matching (stricter)")
    print(f"  3. Combined analysis")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)
    
    return metaproteome_tier1_file, fasta_file

if __name__ == "__main__":
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    metaproteome_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation'
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    create_refined_tier1_for_metaproteome(intensity_analysis_dir, metaproteome_dir, log_dir)