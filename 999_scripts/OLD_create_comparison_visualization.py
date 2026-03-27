import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def create_method_comparison_visualization(output_dir):
    """
    Create comparison visualization with legend below stacked bar, original table position.
    """
    
    print("="*70)
    print("CREATING METHOD COMPARISON VISUALIZATION")
    print("="*70)
    
    # Load all three distribution files
    fuzzy_file = f"{output_dir}/fuzzy_species_distribution.txt"
    exact_file = f"{output_dir}/exact_species_distribution.txt"
    iedb_file = f"{output_dir}/iedb_species_distribution.txt"
    
    print("\nLoading distribution files...")
    fuzzy_df = pd.read_csv(fuzzy_file, sep='\t')
    exact_df = pd.read_csv(exact_file, sep='\t')
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    
    print(f"  Fuzzy: {len(fuzzy_df)} epitopes")
    print(f"  Exact: {len(exact_df)} epitopes")
    print(f"  IEDB: {len(iedb_df)} epitopes")
    
    # Create comparison data
    methods = {
        'Fuzzy K-mer': fuzzy_df,
        'Exact K-mer': exact_df,
        'IEDB BLAST': iedb_df
    }
    
    comparison_data = []
    
    for method_name, df in methods.items():
        total = len(df)
        single = (df['Num_Species'] == 1).sum()
        narrow = ((df['Num_Species'] >= 2) & (df['Num_Species'] <= 3)).sum()
        moderate = ((df['Num_Species'] >= 4) & (df['Num_Species'] <= 6)).sum()
        broad = (df['Num_Species'] >= 7).sum()
        
        comparison_data.append({
            'Method': method_name,
            'Total_Epitopes': total,
            'Single_Species': single,
            'Single_Species_Pct': (single/total)*100,
            'Narrow_2-3': narrow,
            'Narrow_2-3_Pct': (narrow/total)*100,
            'Moderate_4-6': moderate,
            'Moderate_4-6_Pct': (moderate/total)*100,
            'Broad_7plus': broad,
            'Broad_7plus_Pct': (broad/total)*100,
            'Max_Species': int(df['Num_Species'].max()),
            'Median_Species': df['Num_Species'].median()
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison table
    comp_file = f"{output_dir}/method_comparison_species.txt"
    comparison_df.to_csv(comp_file, sep='\t', index=False)
    print(f"\n✓ Saved: {comp_file}")
    
    # Print comparison
    print("\nCOMPARISON TABLE:")
    print(comparison_df.to_string(index=False))
    
    # Create visualization - back to 2x2 layout
    print("\nCreating visualization...")
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    
    colors_methods = {
        'Fuzzy K-mer': '#e78ac3',
        'Exact K-mer': '#8da0cb',
        'IEDB BLAST': '#66c2a5'
    }
    
    # 1. Side-by-side histograms (with log scale for visibility)
    ax1 = axes[0, 0]
    
    max_species = max(df['Num_Species'].max() for df in methods.values())
    x_positions = np.arange(1, max_species + 1)
    width = 0.25
    
    for i, (method_name, df) in enumerate(methods.items()):
        species_counts = df['Num_Species'].value_counts().reindex(x_positions, fill_value=0)
        offset = (i - 1) * width
        ax1.bar(x_positions + offset, species_counts.values, width, 
               label=method_name, color=colors_methods[method_name], edgecolor='black', alpha=0.8)
    
    ax1.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Epitopes (log scale)', fontsize=12, fontweight='bold')
    ax1.set_title('Species Distribution Comparison Across Methods', fontsize=14, fontweight='bold')
    ax1.set_yscale('log')
    ax1.legend(loc='upper right')
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. Stacked bar chart of categories
    ax2 = axes[0, 1]
    
    categories = ['Single_Species_Pct', 'Narrow_2-3_Pct', 'Moderate_4-6_Pct', 'Broad_7plus_Pct']
    category_labels = ['1 species', '2-3 species', '4-6 species', '7+ species']
    colors_cat = ['#fc8d62', '#e78ac3', '#8da0cb', '#66c2a5']
    
    x = np.arange(len(comparison_df))
    bottom = np.zeros(len(comparison_df))
    
    for i, (cat, label) in enumerate(zip(categories, category_labels)):
        values = comparison_df[cat].values
        ax2.bar(x, values, 0.6, label=label, bottom=bottom, color=colors_cat[i], edgecolor='black')
        
        # Add percentage labels in the middle of each segment
        for j, val in enumerate(values):
            if val > 3:  # Only show label if segment is big enough
                ax2.text(j, bottom[j] + val/2, f'{val:.1f}%', 
                        ha='center', va='center', fontsize=9, fontweight='bold')
        
        bottom += values
    
    ax2.set_ylabel('Percentage of Epitopes', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution Categories by Method', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(comparison_df['Method'], fontsize=11)
    # Legend below the plot, horizontal
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4, frameon=True)
    ax2.set_ylim(0, 100)
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Box plot comparison
    ax3 = axes[1, 0]
    
    data_for_box = [df['Num_Species'].values for df in methods.values()]
    labels_for_box = list(methods.keys())
    colors_for_box = [colors_methods[m] for m in methods.keys()]
    
    bp = ax3.boxplot(data_for_box, labels=labels_for_box, patch_artist=True,
                     widths=0.6, showfliers=False)
    
    for patch, color in zip(bp['boxes'], colors_for_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax3.set_ylabel('Number of Species', fontsize=12, fontweight='bold')
    ax3.set_title('Species Count Distribution (Box Plot)', fontsize=14, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # Add median values as text
    for i, (method_name, df) in enumerate(methods.items(), 1):
        median = df['Num_Species'].median()
        ax3.text(i, median, f'{median:.1f}', ha='center', va='bottom', 
                fontweight='bold', fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 4. Summary statistics table (original position, lower title)
    ax4 = axes[1, 1]
    ax4.axis('tight')
    ax4.axis('off')
    
    table_data = []
    for _, row in comparison_df.iterrows():
        table_data.append([
            row['Method'],
            f"{row['Total_Epitopes']:.0f}",
            f"{row['Single_Species_Pct']:.1f}%",
            f"{row['Broad_7plus_Pct']:.1f}%",
            f"{row['Max_Species']:.0f}",
            f"{row['Median_Species']:.1f}"
        ])
    
    table = ax4.table(cellText=table_data,
                     colLabels=['Method', 'Total', '1 Species', '7+ Species', 'Max', 'Median'],
                     cellLoc='center',
                     loc='center',
                     colWidths=[0.25, 0.12, 0.15, 0.15, 0.12, 0.12])
    
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)
    
    # Color header
    for i in range(6):
        table[(0, i)].set_facecolor('#cccccc')
        table[(0, i)].set_text_props(weight='bold')
    
    # Color method names
    for i, method in enumerate(comparison_df['Method'], 1):
        table[(i, 0)].set_facecolor(colors_methods[method])
        table[(i, 0)].set_alpha(0.3)
    
    # Title LOWERED (more padding)
    ax4.set_title('Summary Statistics', fontsize=14, fontweight='bold', pad=50)
    
    plt.suptitle('Method Comparison: Species Distribution Analysis', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    vis_file = f"{output_dir}/method_comparison_visualization.png"
    plt.savefig(vis_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {vis_file}")
    plt.close()
    
    print("\n" + "="*70)
    print("COMPARISON VISUALIZATION COMPLETE")
    print("="*70)

if __name__ == "__main__":
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    create_method_comparison_visualization(output_dir)