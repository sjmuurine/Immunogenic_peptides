import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from datetime import datetime

def visualize_intensity_data_improved(intensity_analysis_dir, output_dir):
    """
    Create improved visualizations with better scaling and zoom.
    """
    
    print("="*80)
    print("ANATOMIA INTENSITY DATA VISUALIZATION (IMPROVED)")
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
    
    # Find intensity columns
    intensity_cols = [col for col in df.columns if col.endswith('_Mean') and 'Intensity' in col]
    print(f"Found intensity columns: {intensity_cols}")
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Main overview with improved scaling
    print("\n1. Creating improved overview plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('ANATOMIA Intensity Analysis for Epitope Clusters', fontsize=16, fontweight='bold')
    
    # Scatter plot: Intensity vs Priority Score (improved)
    if len(intensity_cols) > 0:
        main_intensity = intensity_cols[0]
        
        ax = axes[0, 0]
        for tier in ['TIER1', 'TIER2', 'TIER3']:
            tier_data = anatomia_df[anatomia_df['Tier'] == tier]
            if len(tier_data) > 0:
                ax.scatter(tier_data['Priority_Score'], tier_data[main_intensity], 
                          label=f'{tier} (n={len(tier_data)})', alpha=0.7, s=60)
        
        ax.set_xlabel('Priority Score')
        ax.set_ylabel(main_intensity.replace('_', ' '))
        ax.set_title('Intensity vs Priority Score by Tier')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Improve y-axis scaling
        y_max = anatomia_df[main_intensity].max()
        if y_max > 0:
            ax.set_ylim(-y_max*0.05, y_max*1.1)
    
    # 2. Box plot: Intensity distribution by Tier (fixed)
    ax = axes[0, 1]
    if len(intensity_cols) > 0:
        tier_data_for_plot = []
        tier_labels = []
        
        for tier in ['TIER1', 'TIER2', 'TIER3']:
            tier_subset = anatomia_df[anatomia_df['Tier'] == tier]
            if len(tier_subset) > 0:
                values = tier_subset[main_intensity][tier_subset[main_intensity] > 0]
                if len(values) > 0:
                    tier_data_for_plot.append(values)
                    tier_labels.append(f'{tier}\n(n={len(values)})')
        
        if tier_data_for_plot:
            bp = ax.boxplot(tier_data_for_plot, labels=tier_labels, patch_artist=True)
            # Color the boxes
            colors = ['lightcoral', 'lightblue', 'lightgreen']
            for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                patch.set_facecolor(color)
            
            ax.set_ylabel(main_intensity.replace('_', ' '))
            ax.set_title('Intensity Distribution by Tier')
            ax.grid(True, alpha=0.3)
    
    # 3. Species vs Intensity (improved scaling)
    ax = axes[1, 0]
    if len(intensity_cols) > 0:
        scatter = ax.scatter(anatomia_df['Num_Species'], anatomia_df[main_intensity], 
                   c=anatomia_df['Max_Human_Identity'], cmap='viridis_r', alpha=0.7, s=60)
        ax.set_xlabel('Number of Species')
        ax.set_ylabel(main_intensity.replace('_', ' '))
        ax.set_title('Intensity vs Species Count\n(Color = Human Identity %)')
        
        # Improve scaling
        ax.set_xlim(0.5, anatomia_df['Num_Species'].max() + 0.5)
        y_max = anatomia_df[main_intensity].max()
        if y_max > 0:
            ax.set_ylim(-y_max*0.05, y_max*1.1)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Human Identity %')
        ax.grid(True, alpha=0.3)
    
    # 4. Human Identity vs Intensity (IMPROVED SCALING)
    ax = axes[1, 1]
    if len(intensity_cols) > 0:
        # Find high-intensity epitopes in TIER2/3
        tier1_intensities = anatomia_df[anatomia_df['Tier'] == 'TIER1'][main_intensity]
        if len(tier1_intensities) > 0:
            tier1_75th = np.percentile(tier1_intensities[tier1_intensities > 0], 75)
        else:
            tier1_75th = anatomia_df[main_intensity].quantile(0.75)
        
        high_intensity_lower_tiers = anatomia_df[
            (anatomia_df['Tier'].isin(['TIER2', 'TIER3'])) & 
            (anatomia_df[main_intensity] > tier1_75th)
        ]
        
        # Plot all points
        for tier in ['TIER1', 'TIER2', 'TIER3']:
            tier_data = anatomia_df[anatomia_df['Tier'] == tier]
            if len(tier_data) > 0:
                ax.scatter(tier_data['Max_Human_Identity'], tier_data[main_intensity], 
                          label=tier, alpha=0.6, s=50)
        
        # Highlight high-intensity lower tier candidates
        if len(high_intensity_lower_tiers) > 0:
            ax.scatter(high_intensity_lower_tiers['Max_Human_Identity'], 
                      high_intensity_lower_tiers[main_intensity],
                      s=120, facecolors='none', edgecolors='red', linewidth=2, 
                      label=f'High-intensity TIER2/3 (n={len(high_intensity_lower_tiers)})')
        
        ax.set_xlabel('Max Human Identity %')
        ax.set_ylabel(main_intensity.replace('_', ' '))
        ax.set_title('Human Identity vs Intensity\n(Red circles = high-intensity lower tiers)')
        
        # IMPROVED SCALING - zoom into the relevant range
        identity_min = anatomia_df['Max_Human_Identity'].min()
        identity_max = anatomia_df['Max_Human_Identity'].max()
        
        # Exclude the outlier at 0% for better scaling
        identity_filtered = anatomia_df[anatomia_df['Max_Human_Identity'] > 10]['Max_Human_Identity']
        if len(identity_filtered) > 0:
            ax.set_xlim(identity_filtered.min() - 5, identity_max + 2)
        else:
            ax.set_xlim(identity_min - 5, identity_max + 2)
        
        y_max = anatomia_df[main_intensity].max()
        if y_max > 0:
            ax.set_ylim(-y_max*0.05, y_max*1.1)
        
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    overview_file = f"{output_dir}/intensity_overview_improved.png"
    plt.savefig(overview_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {overview_file}")
    plt.close()
    
    # 2. Detailed intensity comparison (FIXED)
    if len(intensity_cols) > 1:
        print("\n2. Creating detailed intensity comparison...")
        
        n_cols = min(len(intensity_cols), 3)
        n_rows = (len(intensity_cols) + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        
        # Handle single row case
        if n_rows == 1 and n_cols == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        fig.suptitle('Detailed Intensity Analysis Across All Metrics', fontsize=16, fontweight='bold')
        
        for i, col in enumerate(intensity_cols):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # Box plot by tier
            tier_data_for_plot = []
            tier_labels = []
            
            for tier in ['TIER1', 'TIER2', 'TIER3']:
                tier_subset = anatomia_df[anatomia_df['Tier'] == tier]
                if len(tier_subset) > 0:
                    values = tier_subset[col][tier_subset[col] > 0]
                    if len(values) > 0:
                        tier_data_for_plot.append(values)
                        tier_labels.append(f'{tier}')
            
            if tier_data_for_plot:
                bp = ax.boxplot(tier_data_for_plot, labels=tier_labels, patch_artist=True)
                # Color the boxes
                colors = ['lightcoral', 'lightblue', 'lightgreen']
                for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                    patch.set_facecolor(color)
                    
                ax.set_ylabel(col.replace('_', ' '))
                ax.set_title(col.replace('_', ' '))
                ax.grid(True, alpha=0.3)
        
        # Hide empty subplots
        for i in range(len(intensity_cols), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        detailed_file = f"{output_dir}/intensity_detailed_fixed.png"
        plt.savefig(detailed_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {detailed_file}")
        plt.close()
    
    # 3. High-intensity candidates analysis
    print("\n3. Finding high-intensity candidates in lower tiers...")
    
    if len(intensity_cols) > 0:
        main_intensity = intensity_cols[0]
        
        # Calculate 75th percentile threshold
        tier1_intensities = anatomia_df[anatomia_df['Tier'] == 'TIER1'][main_intensity]
        if len(tier1_intensities) > 0:
            threshold = np.percentile(tier1_intensities[tier1_intensities > 0], 75)
        else:
            threshold = anatomia_df[main_intensity].quantile(0.9)
        
        high_intensity_candidates = anatomia_df[
            (anatomia_df['Tier'].isin(['TIER2', 'TIER3'])) & 
            (anatomia_df[main_intensity] > threshold)
        ].sort_values(main_intensity, ascending=False)
        
        print(f"\nHigh-intensity candidates in TIER2/3 (threshold: {threshold:.6f}):")
        print(f"Found {len(high_intensity_candidates)} candidates")
        
        if len(high_intensity_candidates) > 0:
            # Save to file
            candidates_file = f"{output_dir}/high_intensity_lower_tier_candidates.txt"
            high_intensity_candidates[['Cluster_ID', 'Representative_Sequence', 'Tier', 'Num_Species', 
                                     'Max_Human_Identity', 'Priority_Score'] + intensity_cols].to_csv(
                candidates_file, sep='\t', index=False)
            print(f"✓ Saved candidates: {candidates_file}")
            
            # Display top 10
            print(f"\nTop 10 high-intensity lower-tier candidates:")
            print(f"{'Sequence':<35} {'Tier':<6} {'Sp':<3} {'ID%':<6} {'Intensity':<12}")
            print("-" * 80)
            for _, row in high_intensity_candidates.head(10).iterrows():
                seq = row['Representative_Sequence'][:33]
                tier = row['Tier']
                species = row['Num_Species']
                identity = row['Max_Human_Identity']
                intensity = row[main_intensity]
                print(f"{seq:<35} {tier:<6} {species:<3} {identity:<6.1f} {intensity:<12.6f}")
    
    print(f"\n" + "="*80)
    print("IMPROVED VISUALIZATION COMPLETE")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  {overview_file}")
    if len(intensity_cols) > 1:
        print(f"  {detailed_file}")
    if len(high_intensity_candidates) > 0:
        print(f"  {candidates_file}")

if __name__ == "__main__":
    intensity_analysis_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_intensity_analysis'
    output_dir = f"{intensity_analysis_dir}/visualizations"
    
    visualize_intensity_data_improved(intensity_analysis_dir, output_dir)