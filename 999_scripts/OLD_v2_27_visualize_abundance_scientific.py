import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime

def visualize_intensity_data_scientific(intensity_analysis_dir, output_dir):
    """
    Create publication-quality visualizations with proper intensity metrics and colorblind-friendly colors.
    """
    
    print("="*80)
    print("ANATOMIA INTENSITY VISUALIZATION (PUBLICATION QUALITY)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    input_file = f"{intensity_analysis_dir}/ALL_clusters_with_intensity_fixed.txt"
    print(f"\nLoading: {input_file}")
    
    df = pd.read_csv(input_file, sep='\t')
    print(f"Total clusters: {len(df)}")
    
    # Filter to epitopes with ANATOMIA data
    anatomia_df = df[df['Has_ANATOMIA_Data'] == True].copy()
    print(f"Clusters with ANATOMIA data: {len(anatomia_df)}")
    
    if len(anatomia_df) == 0:
        print("No epitopes with ANATOMIA data found!")
        return
    
    # Calculate proper intensity metric: Sum of Protein_Mean_Relative_Intensity
    print("\nCalculating Total Relative Abundance (sum of protein intensities)...")
    
    intensity_file = f"{intensity_analysis_dir}/../001_data/anatomia_protein_intensities.txt"
    intensity_raw_df = pd.read_csv(intensity_file, sep='\t')
    
    # Create protein ID to intensity mapping
    protein_intensity_map = dict(zip(intensity_raw_df['Protein_ID'], 
                                   intensity_raw_df['Protein_Mean_Relative_Intensity']))
    
    def calculate_total_abundance(proteins_str):
        """Calculate sum of protein relative abundances for a cluster."""
        if pd.isna(proteins_str) or proteins_str == '':
            return 0.0
        
        total_abundance = 0.0
        protein_count = 0
        
        for protein_full in proteins_str.split(';'):
            protein_full = protein_full.strip()
            if protein_full.startswith('WP_'):
                protein_id = protein_full.split()[0]  # Extract WP_XXXXX.X
                if protein_id in protein_intensity_map:
                    abundance = protein_intensity_map[protein_id]
                    if pd.notna(abundance) and abundance > 0:
                        total_abundance += abundance
                        protein_count += 1
        
        return total_abundance
    
    # Calculate total abundance for each cluster
    anatomia_df['Total_Relative_Abundance'] = anatomia_df['Source_Proteins_All'].apply(calculate_total_abundance)
    anatomia_df['Has_Abundance_Data'] = anatomia_df['Total_Relative_Abundance'] > 0
    
    # Filter to clusters with abundance data
    abundance_df = anatomia_df[anatomia_df['Has_Abundance_Data']].copy()
    print(f"Clusters with abundance data: {len(abundance_df)}")
    
    if len(abundance_df) == 0:
        print("No clusters with abundance data found!")
        return
    
    # Set up colorblind-friendly style
    plt.style.use('default')
    # Colorblind-friendly palette: blue, orange, purple
    colors = ['#1f77b4', '#ff7f0e', '#9467bd']  # Blue, Orange, Purple
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
    
    # 1. Main overview with proper abundance metrics
    print("\n1. Creating publication-quality overview...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('ANATOMIA Protein Abundance Analysis for Epitope Clusters', fontsize=16, fontweight='bold')
    
    # Scatter plot: Abundance vs Priority Score
    ax = axes[0, 0]
    for i, tier in enumerate(['TIER1', 'TIER2', 'TIER3']):
        tier_data = abundance_df[abundance_df['Tier'] == tier]
        if len(tier_data) > 0:
            ax.scatter(tier_data['Priority_Score'], tier_data['Total_Relative_Abundance'], 
                      label=f'{tier} (n={len(tier_data)})', alpha=0.7, s=60, color=colors[i])
    
    ax.set_xlabel('Priority Score')
    ax.set_ylabel('Total Relative Abundance')
    ax.set_title('Protein Abundance vs Priority Score by Tier')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Improve y-axis scaling
    y_max = abundance_df['Total_Relative_Abundance'].max()
    if y_max > 0:
        ax.set_ylim(-y_max*0.05, y_max*1.1)
    
    # 2. Box plot: Abundance distribution by Tier
    ax = axes[0, 1]
    tier_data_for_plot = []
    tier_labels = []
    
    for tier in ['TIER1', 'TIER2', 'TIER3']:
        tier_subset = abundance_df[abundance_df['Tier'] == tier]
        if len(tier_subset) > 0:
            values = tier_subset['Total_Relative_Abundance']
            if len(values) > 0:
                tier_data_for_plot.append(values)
                tier_labels.append(f'{tier}\n(n={len(values)})')
    
    if tier_data_for_plot:
        bp = ax.boxplot(tier_data_for_plot, labels=tier_labels, patch_artist=True)
        # Color the boxes with colorblind-friendly colors
        for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Total Relative Abundance')
        ax.set_title('Abundance Distribution by Tier')
        ax.grid(True, alpha=0.3)
    
    # 3. Species vs Abundance (colorblind-safe colormap)
    ax = axes[1, 0]
    scatter = ax.scatter(abundance_df['Num_Species'], abundance_df['Total_Relative_Abundance'], 
               c=abundance_df['Max_Human_Identity'], cmap='plasma', alpha=0.7, s=60)  # plasma is colorblind-safe
    ax.set_xlabel('Number of Species')
    ax.set_ylabel('Total Relative Abundance')
    ax.set_title('Abundance vs Species Count\n(Color = Human Identity %)')
    
    # Improve scaling
    ax.set_xlim(0.5, abundance_df['Num_Species'].max() + 0.5)
    y_max = abundance_df['Total_Relative_Abundance'].max()
    if y_max > 0:
        ax.set_ylim(-y_max*0.05, y_max*1.1)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Human Identity %')
    ax.grid(True, alpha=0.3)
    
    # 4. Human Identity vs Abundance (improved scaling + colorblind colors)
    ax = axes[1, 1]
    
    # Find high-abundance epitopes in TIER2/3
    tier1_abundances = abundance_df[abundance_df['Tier'] == 'TIER1']['Total_Relative_Abundance']
    if len(tier1_abundances) > 0:
        threshold = np.percentile(tier1_abundances, 75)
    else:
        threshold = abundance_df['Total_Relative_Abundance'].quantile(0.75)
    
    high_abundance_lower_tiers = abundance_df[
        (abundance_df['Tier'].isin(['TIER2', 'TIER3'])) & 
        (abundance_df['Total_Relative_Abundance'] > threshold)
    ]
    
    # Plot all points with colorblind-friendly colors
    for i, tier in enumerate(['TIER1', 'TIER2', 'TIER3']):
        tier_data = abundance_df[abundance_df['Tier'] == tier]
        if len(tier_data) > 0:
            ax.scatter(tier_data['Max_Human_Identity'], tier_data['Total_Relative_Abundance'], 
                      label=tier, alpha=0.6, s=50, color=colors[i])
    
    # Highlight high-abundance lower tier candidates with distinct marker
    if len(high_abundance_lower_tiers) > 0:
        ax.scatter(high_abundance_lower_tiers['Max_Human_Identity'], 
                  high_abundance_lower_tiers['Total_Relative_Abundance'],
                  s=120, facecolors='none', edgecolors='black', linewidth=3, marker='^',
                  label=f'High-abundance TIER2/3 (n={len(high_abundance_lower_tiers)})')
    
    ax.set_xlabel('Max Human Identity %')
    ax.set_ylabel('Total Relative Abundance')
    ax.set_title('Human Identity vs Abundance\n(Triangles = high-abundance lower tiers)')
    
    # Improved scaling - zoom into relevant range
    identity_filtered = abundance_df[abundance_df['Max_Human_Identity'] > 10]['Max_Human_Identity']
    if len(identity_filtered) > 0:
        ax.set_xlim(identity_filtered.min() - 5, abundance_df['Max_Human_Identity'].max() + 2)
    else:
        ax.set_xlim(abundance_df['Max_Human_Identity'].min() - 5, abundance_df['Max_Human_Identity'].max() + 2)
    
    y_max = abundance_df['Total_Relative_Abundance'].max()
    if y_max > 0:
        ax.set_ylim(-y_max*0.05, y_max*1.1)
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    overview_file = f"{output_dir}/abundance_overview_publication.png"
    plt.savefig(overview_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {overview_file}")
    plt.close()
    
    # 2. High-abundance candidates analysis and table
    print("\n2. Identifying high-abundance candidates in lower tiers...")
    
    high_abundance_candidates = abundance_df[
        (abundance_df['Tier'].isin(['TIER2', 'TIER3'])) & 
        (abundance_df['Total_Relative_Abundance'] > threshold)
    ].sort_values('Total_Relative_Abundance', ascending=False)
    
    print(f"\nHigh-abundance candidates in TIER2/3 (threshold: {threshold:.4f}):")
    print(f"Found {len(high_abundance_candidates)} candidates")
    
    if len(high_abundance_candidates) > 0:
        # Save to file
        candidates_file = f"{output_dir}/high_abundance_lower_tier_candidates_scientific.txt"
        output_cols = ['Cluster_ID', 'Representative_Sequence', 'Tier', 'Num_Species', 
                      'Max_Human_Identity', 'Priority_Score', 'Total_Relative_Abundance', 
                      'Num_Matched_ANATOMIA', 'Species_List']
        high_abundance_candidates[output_cols].to_csv(candidates_file, sep='\t', index=False)
        print(f"✓ Saved candidates: {candidates_file}")
        
        # Display top 15
        print(f"\nTop 15 high-abundance lower-tier candidates:")
        print(f"{'Sequence':<35} {'Tier':<6} {'Sp':<3} {'ID%':<6} {'Abundance':<12}")
        print("-" * 80)
        for _, row in high_abundance_candidates.head(15).iterrows():
            seq = row['Representative_Sequence'][:33]
            tier = row['Tier']
            species = row['Num_Species']
            identity = row['Max_Human_Identity']
            abundance = row['Total_Relative_Abundance']
            print(f"{seq:<35} {tier:<6} {species:<3} {identity:<6.1f} {abundance:<12.4f}")
    
    # 3. Summary statistics
    print(f"\n" + "="*80)
    print("ABUNDANCE ANALYSIS SUMMARY")
    print("="*80)
    
    print(f"\nIntensity metric used: Total Relative Abundance")
    print(f"  = Sum of Protein_Mean_Relative_Intensity for all proteins in cluster")
    print(f"  Publication description: 'Total relative protein abundance'")
    
    print(f"\nData overview:")
    print(f"  Total clusters: {len(df)}")
    print(f"  With ANATOMIA data: {len(anatomia_df)} ({len(anatomia_df)/len(df)*100:.1f}%)")
    print(f"  With abundance data: {len(abundance_df)} ({len(abundance_df)/len(df)*100:.1f}%)")
    
    for tier in ['TIER1', 'TIER2', 'TIER3']:
        tier_abundance = abundance_df[abundance_df['Tier'] == tier]
        if len(tier_abundance) > 0:
            mean_abundance = tier_abundance['Total_Relative_Abundance'].mean()
            print(f"  {tier} mean abundance: {mean_abundance:.4f}")
    
    print(f"\nColorblind-friendly features:")
    print(f"  ✓ Blue/Orange/Purple color scheme")
    print(f"  ✓ Plasma colormap (colorblind-safe)")
    print(f"  ✓ Distinct markers (triangles vs circles)")
    
    print(f"\nOutput files:")
    print(f"  Main visualization: {overview_file}")
    if len(high_abundance_candidates) > 0:
        print(f"  High-abundance candidates: {candidates_file}")
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)

if __name__ == "__main__":
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    output_dir = f"{intensity_analysis_dir}/visualizations"
    
    visualize_intensity_data_scientific(intensity_analysis_dir, output_dir)