import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def parse_species_from_header(header_str):
    """Parse species from ANATOMIA or LITERATURE header."""
    if pd.isna(header_str):
        return None
    
    header_str = str(header_str)
    
    # LITERATURE format
    if 'LITERATURE' in header_str:
        species_map = {
            'fragilis': 'B_fragilis', 'nucleatum': 'F_nucleatum', 'muciniphila': 'A_muciniphila',
            'prausnitzii': 'F_prausnitzii', 'reuteri': 'L_reuteri', 'cloacae': 'E_cloacae',
            'faecalis': 'E_faecalis', 'salivarius': 'S_salivarius', 'longum': 'B_longum',
            'plantarum': 'L_plantarum', 'casei': 'L_casei', 'sanguinis': 'T_sanguinis', 
            'magna': 'V_magna', 'russellii': 'P_russellii'
        }
        for key, value in species_map.items():
            if key in header_str.lower():
                return value
        return 'LITERATURE'
    
    # ANATOMIA format: "ProteinID Species_OG_Conservation_Localization"
    if ' ' not in header_str:
        return None
    
    parts = header_str.split(' ', 1)[1].split('_')
    
    # Find OG position
    og_idx = None
    for i, part in enumerate(parts):
        if part.startswith('OG') or part == 'unique':
            og_idx = i
            break
    
    if og_idx and og_idx > 0:
        return '_'.join(parts[:og_idx])
    
    return None

def analyze_exact_kmer(output_dir):
    """Analyze exact k-mer species distribution."""
    
    print("\n" + "="*70)
    print("ANALYZING EXACT K-MER")
    print("="*70)
    
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    exact_df = pd.read_csv(exact_file, sep='\t')
    
    print(f"Total exact k-mer motifs: {len(exact_df)}")
    
    # Parse species from Sequence_IDs (contains multiple protein headers separated by "; ")
    exact_results = []
    
    for _, row in exact_df.iterrows():
        motif = row['Motif']
        seq_ids = str(row['Sequence_IDs'])
        
        # Split by "; " to get individual protein headers
        headers = seq_ids.split('; ')
        
        species_set = set()
        for header in headers:
            species = parse_species_from_header(header)
            if species and species != 'LITERATURE':
                species_set.add(species)
        
        exact_results.append({
            'Epitope_Sequence': motif,
            'Species_List': sorted(list(species_set)),
            'Total_Hits': row['UniqueSeqCount']
        })
    
    exact_summary = pd.DataFrame(exact_results)
    exact_summary['Num_Species'] = exact_summary['Species_List'].apply(len)
    exact_summary['Species_Names'] = exact_summary['Species_List'].apply(lambda x: '; '.join(x))
    exact_summary = exact_summary.sort_values('Num_Species', ascending=True)
    
    # Save
    exact_output = f"{output_dir}/exact_species_distribution.txt"
    exact_summary[['Epitope_Sequence', 'Num_Species', 'Total_Hits', 'Species_Names']].to_csv(
        exact_output, sep='\t', index=False
    )
    print(f"✓ Saved: {exact_output}")
    
    return exact_summary

def analyze_iedb_blast(output_dir):
    """Analyze IEDB BLAST species distribution."""
    
    print("\n" + "="*70)
    print("ANALYZING IEDB BLAST")
    print("="*70)
    
    iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    
    print(f"Total IEDB BLAST records: {len(iedb_df)}")
    
    # Parse species from Protein_Full_Header
    iedb_df['Species'] = iedb_df['Protein_Full_Header'].apply(parse_species_from_header)
    
    # Filter to ANATOMIA only
    iedb_anatomia = iedb_df[iedb_df['Data_Source'] == 'ANATOMIA'].copy()
    
    print(f"ANATOMIA records: {len(iedb_anatomia)}")
    print(f"Species parsed: {iedb_anatomia['Species'].notna().sum()}")
    
    # Group by epitope sequence
    iedb_summary = iedb_anatomia.groupby('Epitope_Sequence').agg({
        'Species': lambda x: sorted(list(x.dropna().unique())),
        'Protein_Full_Header': 'count'
    }).reset_index()
    
    iedb_summary.columns = ['Epitope_Sequence', 'Species_List', 'Total_Hits']
    iedb_summary['Num_Species'] = iedb_summary['Species_List'].apply(len)
    iedb_summary['Species_Names'] = iedb_summary['Species_List'].apply(lambda x: '; '.join(x))
    iedb_summary = iedb_summary.sort_values('Num_Species', ascending=True)
    
    # Save
    iedb_output = f"{output_dir}/iedb_species_distribution.txt"
    iedb_summary[['Epitope_Sequence', 'Num_Species', 'Total_Hits', 'Species_Names']].to_csv(
        iedb_output, sep='\t', index=False
    )
    print(f"✓ Saved: {iedb_output}")
    
    return iedb_summary

def create_method_comparison(fuzzy_df, exact_df, iedb_df, output_dir):
    """Create comparison across all three methods."""
    
    print("\n" + "="*70)
    print("CREATING METHOD COMPARISON")
    print("="*70)
    
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
            'Max_Species': df['Num_Species'].max(),
            'Median_Species': df['Num_Species'].median()
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save
    comp_file = f"{output_dir}/method_comparison_species.txt"
    comparison_df.to_csv(comp_file, sep='\t', index=False)
    print(f"✓ Saved: {comp_file}")
    
    # Print summary
    print("\nCOMPARISON SUMMARY:")
    print(comparison_df.to_string(index=False))
    
    return comparison_df

def visualize_all_methods(fuzzy_df, exact_df, iedb_df, comparison_df, output_dir):
    """Create comparative visualizations."""
    
    print("\n" + "="*70)
    print("CREATING COMPARATIVE VISUALIZATIONS")
    print("="*70)
    
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    
    methods = {
        'Fuzzy K-mer': (fuzzy_df, '#e78ac3'),
        'Exact K-mer': (exact_df, '#8da0cb'),
        'IEDB BLAST': (iedb_df, '#66c2a5')
    }
    
    # 1. Side-by-side histograms
    ax1 = axes[0, 0]
    
    max_species = max(df['Num_Species'].max() for df, _ in methods.values())
    x_positions = np.arange(1, max_species + 1)
    width = 0.25
    
    for i, (method_name, (df, color)) in enumerate(methods.items()):
        species_counts = df['Num_Species'].value_counts().reindex(x_positions, fill_value=0)
        offset = (i - 1) * width
        ax1.bar(x_positions + offset, species_counts.values, width, 
               label=method_name, color=color, edgecolor='black', alpha=0.8)
    
    ax1.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Epitopes', fontsize=12, fontweight='bold')
    ax1.set_title('Species Distribution Comparison Across Methods', fontsize=14, fontweight='bold')
    ax1.legend()
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
        bottom += values
    
    ax2.set_ylabel('Percentage of Epitopes', fontsize=12, fontweight='bold')
    ax2.set_title('Distribution Categories by Method', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(comparison_df['Method'])
    ax2.legend(loc='upper right')
    ax2.grid(axis='y', alpha=0.3)
    
    # 3. Box plot comparison
    ax3 = axes[1, 0]
    
    data_for_box = [df['Num_Species'].values for df, _ in methods.values()]
    labels_for_box = list(methods.keys())
    colors_for_box = [color for _, color in methods.values()]
    
    bp = ax3.boxplot(data_for_box, labels=labels_for_box, patch_artist=True,
                     widths=0.6, showfliers=False)
    
    for patch, color in zip(bp['boxes'], colors_for_box):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax3.set_ylabel('Number of Species', fontsize=12, fontweight='bold')
    ax3.set_title('Species Count Distribution (Box Plot)', fontsize=14, fontweight='bold')
    ax3.grid(axis='y', alpha=0.3)
    
    # Add median values as text
    for i, (method_name, (df, _)) in enumerate(methods.items(), 1):
        median = df['Num_Species'].median()
        ax3.text(i, median, f'{median:.1f}', ha='center', va='bottom', 
                fontweight='bold', fontsize=10)
    
    # 4. Summary statistics table
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
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Color header
    for i in range(6):
        table[(0, i)].set_facecolor('#cccccc')
        table[(0, i)].set_text_props(weight='bold')
    
    ax4.set_title('Summary Statistics', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    vis_file = f"{output_dir}/method_comparison_visualization.png"
    plt.savefig(vis_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {vis_file}")
    plt.close()

def main():
    """Main analysis function."""
    
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    
    print("="*70)
    print("SPECIES DISTRIBUTION ANALYSIS: ALL METHODS")
    print("="*70)
    
    # Load fuzzy results (already created)
    fuzzy_file = f"{output_dir}/fuzzy_species_distribution.txt"
    print(f"\nLoading existing fuzzy results: {fuzzy_file}")
    fuzzy_df = pd.read_csv(fuzzy_file, sep='\t')
    
    # Analyze exact k-mer
    exact_df = analyze_exact_kmer(output_dir)
    
    # Analyze IEDB BLAST
    iedb_df = analyze_iedb_blast(output_dir)
    
    # Create comparison
    comparison_df = create_method_comparison(fuzzy_df, exact_df, iedb_df, output_dir)
    
    # Create visualizations
    visualize_all_methods(fuzzy_df, exact_df, iedb_df, comparison_df, output_dir)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nFiles created in: {output_dir}")
    print("  - exact_species_distribution.txt")
    print("  - iedb_species_distribution.txt")
    print("  - method_comparison_species.txt")
    print("  - method_comparison_visualization.png")

if __name__ == "__main__":
    main()