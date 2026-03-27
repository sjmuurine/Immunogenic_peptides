import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def visualize_fuzzy_species_distribution(output_dir):
    """
    Visualize and summarize fuzzy epitope species distribution.
    """
    
    print("="*70)
    print("VISUALIZING FUZZY SPECIES DISTRIBUTION")
    print("="*70)
    
    # Load the results we just created
    fuzzy_file = f"{output_dir}/fuzzy_species_distribution.txt"
    
    print(f"\nLoading: {fuzzy_file}")
    fuzzy_df = pd.read_csv(fuzzy_file, sep='\t')
    
    print(f"Total epitopes: {len(fuzzy_df)}")
    
    # Create summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    
    total_epitopes = len(fuzzy_df)
    
    summary_stats = {
        'Category': [],
        'Count': [],
        'Percentage': []
    }
    
    categories = [
        ('1 species (very narrow)', 1, 1),
        ('2-3 species (narrow)', 2, 3),
        ('4-6 species (moderate)', 4, 6),
        ('7+ species (broad)', 7, 100)
    ]
    
    for label, min_sp, max_sp in categories:
        count = ((fuzzy_df['Num_Species'] >= min_sp) & (fuzzy_df['Num_Species'] <= max_sp)).sum()
        pct = (count / total_epitopes) * 100
        summary_stats['Category'].append(label)
        summary_stats['Count'].append(count)
        summary_stats['Percentage'].append(pct)
        print(f"{label:30} {count:6} epitopes ({pct:5.1f}%)")
    
    summary_df = pd.DataFrame(summary_stats)
    summary_file = f"{output_dir}/fuzzy_species_summary.txt"
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"\n✓ Saved summary: {summary_file}")
    
    # Extract epitopes with 4+ species
    broad_epitopes = fuzzy_df[fuzzy_df['Num_Species'] >= 4].copy()
    broad_epitopes = broad_epitopes.sort_values('Num_Species', ascending=False)
    
    print(f"\n" + "="*70)
    print(f"EPITOPES WITH 4+ SPECIES ({len(broad_epitopes)} total)")
    print("="*70)
    
    for _, row in broad_epitopes.iterrows():
        print(f"{row['Epitope_Sequence']:22} - {row['Num_Species']:2} species ({row['Total_Hits']:3} hits): {row['Species_Names']}")
    
    broad_file = f"{output_dir}/fuzzy_broad_epitopes.txt"
    broad_epitopes.to_csv(broad_file, sep='\t', index=False)
    print(f"\n✓ Saved broad epitopes: {broad_file}")
    
    # Create visualizations
    print("\n" + "="*70)
    print("CREATING VISUALIZATIONS")
    print("="*70)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Histogram of species counts
    ax1 = axes[0, 0]
    species_counts = fuzzy_df['Num_Species'].value_counts().sort_index()
    ax1.bar(species_counts.index, species_counts.values, color='#66c2a5', edgecolor='black')
    ax1.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Epitopes', fontsize=12, fontweight='bold')
    ax1.set_title('Fuzzy K-mer: Distribution of Species per Epitope', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add text annotations
    for idx, val in zip(species_counts.index, species_counts.values):
        ax1.text(idx, val, str(val), ha='center', va='bottom', fontsize=9)
    
    # 2. Pie chart of categories
    ax2 = axes[0, 1]
    colors = ['#fc8d62', '#e78ac3', '#8da0cb', '#66c2a5']
    ax2.pie(summary_stats['Count'], labels=summary_stats['Category'], autopct='%1.1f%%',
            colors=colors, startangle=90, textprops={'fontsize': 10})
    ax2.set_title('Fuzzy K-mer: Species Distribution Categories', fontsize=14, fontweight='bold')
    
    # 3. Cumulative distribution
    ax3 = axes[1, 0]
    sorted_counts = fuzzy_df['Num_Species'].sort_values()
    cumulative = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts) * 100
    ax3.plot(sorted_counts, cumulative, linewidth=2, color='#66c2a5')
    ax3.fill_between(sorted_counts, cumulative, alpha=0.3, color='#66c2a5')
    ax3.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Cumulative Percentage (%)', fontsize=12, fontweight='bold')
    ax3.set_title('Fuzzy K-mer: Cumulative Distribution', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% line')
    ax3.axhline(y=90, color='orange', linestyle='--', alpha=0.5, label='90% line')
    ax3.legend()
    
    # 4. Top species in broad epitopes (4+ species)
    ax4 = axes[1, 1]
    if len(broad_epitopes) > 0:
        # Count which species appear most in broad epitopes
        all_species = []
        for species_list in broad_epitopes['Species_Names']:
            all_species.extend(species_list.split('; '))
        
        species_freq = pd.Series(all_species).value_counts().head(10)
        ax4.barh(range(len(species_freq)), species_freq.values, color='#8da0cb', edgecolor='black')
        ax4.set_yticks(range(len(species_freq)))
        ax4.set_yticklabels(species_freq.index, fontsize=10)
        ax4.set_xlabel('Number of Epitopes (4+ species)', fontsize=12, fontweight='bold')
        ax4.set_title('Top Species in Broad-Distribution Epitopes', fontsize=14, fontweight='bold')
        ax4.grid(axis='x', alpha=0.3)
        ax4.invert_yaxis()
        
        for i, val in enumerate(species_freq.values):
            ax4.text(val, i, f' {val}', va='center', fontsize=9)
    else:
        ax4.text(0.5, 0.5, 'No epitopes with 4+ species', 
                ha='center', va='center', fontsize=12)
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
    
    plt.tight_layout()
    
    vis_file = f"{output_dir}/fuzzy_species_distribution.png"
    plt.savefig(vis_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved visualization: {vis_file}")
    plt.close()
    
    # Create comparison table with other methods
    print("\n" + "="*70)
    print("CREATING METHOD COMPARISON")
    print("="*70)
    
    comparison_data = {
        'Method': ['Fuzzy K-mer'],
        'Total_Epitopes': [len(fuzzy_df)],
        'Single_Species': [summary_stats['Count'][0]],
        'Single_Species_Pct': [summary_stats['Percentage'][0]],
        'Narrow_2-3': [summary_stats['Count'][1]],
        'Narrow_2-3_Pct': [summary_stats['Percentage'][1]],
        'Moderate_4-6': [summary_stats['Count'][2]],
        'Moderate_4-6_Pct': [summary_stats['Percentage'][2]],
        'Broad_7plus': [summary_stats['Count'][3]],
        'Broad_7plus_Pct': [summary_stats['Percentage'][3]],
        'Max_Species': [fuzzy_df['Num_Species'].max()]
    }
    
    comparison_df = pd.DataFrame(comparison_data)
    comparison_file = f"{output_dir}/fuzzy_method_comparison.txt"
    comparison_df.to_csv(comparison_file, sep='\t', index=False)
    print(f"✓ Saved comparison: {comparison_file}")
    
    print("\n" + "="*70)
    print("VISUALIZATION COMPLETE")
    print("="*70)
    print(f"\nKey findings:")
    print(f"  - {summary_stats['Percentage'][0]:.1f}% of fuzzy epitopes are species-specific (1 species)")
    print(f"  - {summary_stats['Percentage'][3]:.1f}% have broad distribution (7+ species)")
    print(f"  - Maximum species count: {fuzzy_df['Num_Species'].max()}")
    print(f"  - Median species count: {fuzzy_df['Num_Species'].median():.0f}")

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
visualize_fuzzy_species_distribution(output_dir)