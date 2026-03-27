import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime

def visualize_abundance_data_improved_scales(intensity_analysis_dir, output_dir):
    """
    Create visualizations with properly adjusted scales for better visibility.
    """
    
    print("="*80)
    print("ANATOMIA ABUNDANCE VISUALIZATION (IMPROVED SCALES)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and process data (same as before)
    input_file = f"{intensity_analysis_dir}/ALL_clusters_with_intensity_fixed.txt"
    df = pd.read_csv(input_file, sep='\t')
    anatomia_df = df[df['Has_ANATOMIA_Data'] == True].copy()
    
    # Calculate total abundance
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
    
    # Set up colorblind-friendly colors
    colors = ['#1f77b4', '#ff7f0e', '#9467bd']  # Blue, Orange, Purple
    
    # Create improved visualization
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('ANATOMIA Protein Abundance Analysis (Improved Scales)', fontsize=16, fontweight='bold')
    
    # 1. Scatter plot: Abundance vs Priority Score (unchanged - this one was fine)
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
    
    # 2. Box plot: IMPROVED Y-AXIS SCALING
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
        for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Total Relative Abundance')
        ax.set_title('Abundance Distribution by Tier')
        ax.grid(True, alpha=0.3)
        
        # IMPROVED: Set y-axis to focus on the main data range (exclude extreme outliers)
        all_values = [val for sublist in tier_data_for_plot for val in sublist]
        q99 = np.percentile(all_values, 99)
        ax.set_ylim(0, q99 * 1.1)  # Show up to 99th percentile + 10%
    
    # 3. Species vs Abundance: FIXED COLOR SCALE
    ax = axes[1, 0]
    
    # Calculate proper color range (exclude the 0% outlier)
    identity_data = abundance_df['Max_Human_Identity']
    identity_filtered = identity_data[identity_data > 10]  # Exclude outliers
    vmin = identity_filtered.min() - 2
    vmax = identity_filtered.max() + 2
    
    scatter = ax.scatter(abundance_df['Num_Species'], abundance_df['Total_Relative_Abundance'], 
               c=abundance_df['Max_Human_Identity'], cmap='plasma', alpha=0.7, s=80,
               vmin=vmin, vmax=vmax)  # FIXED: Proper color range
    
    ax.set_xlabel('Number of Species')
    ax.set_ylabel('Total Relative Abundance')
    ax.set_title('Abundance vs Species Count\n(Color = Human Identity %, focused range)')
    
    ax.set_xlim(0.5, abundance_df['Num_Species'].max() + 0.5)
    # Improved y-axis (same as box plot)
    all_abundance = abundance_df['Total_Relative_Abundance']
    q99_abundance = np.percentile(all_abundance, 99)
    ax.set_ylim(0, q99_abundance * 1.1)
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Human Identity %')
    ax.grid(True, alpha=0.3)
    
    # 4. Human Identity vs Abundance: BETTER SEPARATION OF INTERESTING CANDIDATES
    ax = axes[1, 1]
    
    # Find high-abundance candidates in lower tiers
    tier1_abundances = abundance_df[abundance_df['Tier'] == 'TIER1']['Total_Relative_Abundance']
    if len(tier1_abundances) > 0:
        threshold = np.percentile(tier1_abundances, 75)
    else:
        threshold = abundance_df['Total_Relative_Abundance'].quantile(0.75)
    
    high_abundance_lower_tiers = abundance_df[
        (abundance_df['Tier'].isin(['TIER2', 'TIER3'])) & 
        (abundance_df['Total_Relative_Abundance'] > threshold)
    ]
    
    # Plot regular points smaller and more transparent
    for i, tier in enumerate(['TIER1', 'TIER2', 'TIER3']):
        tier_data = abundance_df[abundance_df['Tier'] == tier]
        if len(tier_data) > 0:
            ax.scatter(tier_data['Max_Human_Identity'], tier_data['Total_Relative_Abundance'], 
                      label=tier, alpha=0.4, s=30, color=colors[i])  # SMALLER, more transparent
    
    # Highlight interesting candidates with MUCH larger, more visible markers
    if len(high_abundance_lower_tiers) > 0:
        ax.scatter(high_abundance_lower_tiers['Max_Human_Identity'], 
                  high_abundance_lower_tiers['Total_Relative_Abundance'],
                  s=200, facecolors='none', edgecolors='black', linewidth=4, marker='s',  # MUCH larger squares
                  label=f'High-abundance TIER2/3 (n={len(high_abundance_lower_tiers)})')
        
        # Also add colored fill to make them even more visible
        for _, row in high_abundance_lower_tiers.iterrows():
            tier_color = colors[1] if row['Tier'] == 'TIER2' else colors[2]
            ax.scatter(row['Max_Human_Identity'], row['Total_Relative_Abundance'],
                      s=150, facecolors=tier_color, alpha=0.8, marker='s', edgecolors='black', linewidth=2)
    
    ax.set_xlabel('Max Human Identity %')
    ax.set_ylabel('Total Relative Abundance')
    ax.set_title('Human Identity vs Abundance\n(Large squares = high-abundance lower tiers)')
    
    # Focus on relevant ranges
    ax.set_xlim(vmin, vmax)  # Same as color scale
    ax.set_ylim(0, q99_abundance * 1.1)  # Same as other plots
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    overview_file = f"{output_dir}/abundance_overview_improved_scales.png"
    plt.savefig(overview_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {overview_file}")
    plt.close()
    
    # Print info about the interesting candidates
    if len(high_abundance_lower_tiers) > 0:
        print(f"\nHigh-abundance TIER2/3 candidates (marked with large squares):")
        print(f"{'Sequence':<35} {'Tier':<6} {'Sp':<3} {'ID%':<6} {'Abundance':<12}")
        print("-" * 80)
        for _, row in high_abundance_lower_tiers.sort_values('Total_Relative_Abundance', ascending=False).iterrows():
            seq = row['Representative_Sequence'][:33]
            tier = row['Tier']
            species = row['Num_Species']
            identity = row['Max_Human_Identity']
            abundance = row['Total_Relative_Abundance']
            print(f"{seq:<35} {tier:<6} {species:<3} {identity:<6.1f} {abundance:<12.4f}")
    
    print(f"\nScale improvements:")
    print(f"  ✓ Color range: {vmin:.1f}-{vmax:.1f}% (excluding outliers)")
    print(f"  ✓ Y-axis max: {q99_abundance:.4f} (99th percentile)")
    print(f"  ✓ Large squares highlight high-abundance TIER2/3 candidates")
    print(f"  ✓ Reduced transparency on regular points for better contrast")

if __name__ == "__main__":
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    output_dir = f"{intensity_analysis_dir}/visualizations"
    
    visualize_abundance_data_improved_scales(intensity_analysis_dir, output_dir)