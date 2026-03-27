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
        return None  # Changed from 'LITERATURE' to None to exclude from species counts
    
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

def visualize_single_method(method_name, species_df, output_dir, color='#66c2a5'):
    """
    Create detailed 4-panel visualization for a single method (like fuzzy).
    """
    
    print(f"\nCreating visualization for {method_name}...")
    
    # Calculate statistics
    total_epitopes = len(species_df)
    
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
        count = ((species_df['Num_Species'] >= min_sp) & (species_df['Num_Species'] <= max_sp)).sum()
        pct = (count / total_epitopes) * 100
        summary_stats['Category'].append(label)
        summary_stats['Count'].append(count)
        summary_stats['Percentage'].append(pct)
    
    # Extract broad epitopes (4+ species)
    broad_epitopes = species_df[species_df['Num_Species'] >= 4].copy()
    broad_epitopes = broad_epitopes.sort_values('Num_Species', ascending=False)
    
    # Save broad epitopes
    broad_file = f"{output_dir}/{method_name.lower().replace(' ', '_')}_broad_epitopes.txt"
    broad_epitopes.to_csv(broad_file, sep='\t', index=False)
    
    # Create visualizations
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Histogram of species counts
    ax1 = axes[0, 0]
    species_counts = species_df['Num_Species'].value_counts().sort_index()
    ax1.bar(species_counts.index, species_counts.values, color=color, edgecolor='black')
    ax1.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Epitopes', fontsize=12, fontweight='bold')
    ax1.set_title(f'{method_name}: Distribution of Species per Epitope', fontsize=14, fontweight='bold')
    ax1.grid(axis='y', alpha=0.3)
    
    # Add text annotations
    for idx, val in zip(species_counts.index, species_counts.values):
        ax1.text(idx, val, str(val), ha='center', va='bottom', fontsize=9)
    
    # 2. Pie chart of categories
    ax2 = axes[0, 1]
    colors_pie = ['#fc8d62', '#e78ac3', '#8da0cb', '#66c2a5']
    ax2.pie(summary_stats['Count'], labels=summary_stats['Category'], autopct='%1.1f%%',
            colors=colors_pie, startangle=90, textprops={'fontsize': 10})
    ax2.set_title(f'{method_name}: Species Distribution Categories', fontsize=14, fontweight='bold')
    
    # 3. Cumulative distribution
    ax3 = axes[1, 0]
    sorted_counts = species_df['Num_Species'].sort_values()
    cumulative = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts) * 100
    ax3.plot(sorted_counts, cumulative, linewidth=2, color=color)
    ax3.fill_between(sorted_counts, cumulative, alpha=0.3, color=color)
    ax3.set_xlabel('Number of Species', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Cumulative Percentage (%)', fontsize=12, fontweight='bold')
    ax3.set_title(f'{method_name}: Cumulative Distribution', fontsize=14, fontweight='bold')
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
            if pd.notna(species_list):
                all_species.extend(species_list.split('; '))
        
        if len(all_species) > 0:
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
            ax4.text(0.5, 0.5, 'No species data available', 
                    ha='center', va='center', fontsize=12)
            ax4.set_xlim(0, 1)
            ax4.set_ylim(0, 1)
            ax4.axis('off')
    else:
        ax4.text(0.5, 0.5, 'No epitopes with 4+ species', 
                ha='center', va='center', fontsize=12)
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
    
    plt.tight_layout()
    
    vis_file = f"{output_dir}/{method_name.lower().replace(' ', '_')}_species_visualization.png"
    plt.savefig(vis_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {vis_file}")
    plt.close()
    
    # Print statistics
    print(f"\n{method_name} STATISTICS:")
    print(f"  Total epitopes: {total_epitopes}")
    print(f"  1 species: {summary_stats['Count'][0]} ({summary_stats['Percentage'][0]:.1f}%)")
    print(f"  2-3 species: {summary_stats['Count'][1]} ({summary_stats['Percentage'][1]:.1f}%)")
    print(f"  4-6 species: {summary_stats['Count'][2]} ({summary_stats['Percentage'][2]:.1f}%)")
    print(f"  7+ species: {summary_stats['Count'][3]} ({summary_stats['Percentage'][3]:.1f}%)")
    print(f"  Max species: {species_df['Num_Species'].max()}")
    print(f"  Median species: {species_df['Num_Species'].median():.1f}")

def analyze_exact_kmer(output_dir):
    """Analyze exact k-mer with proper species parsing."""
    
    print("\n" + "="*70)
    print("ANALYZING EXACT K-MER (DETAILED)")
    print("="*70)
    
    exact_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/ANTTIsResults_copies/combined_exact.txt'
    exact_df = pd.read_csv(exact_file, sep='\t')
    
    print(f"Total exact k-mer motifs: {len(exact_df)}")
    
    # Parse species from Sequence_IDs
    exact_results = []
    
    for idx, row in exact_df.iterrows():
        if idx % 100 == 0:
            print(f"  Processing motif {idx+1}/{len(exact_df)}...", end='\r')
        
        motif = row['Motif']
        seq_ids = str(row['Sequence_IDs'])
        
        # Split by "; " to get individual protein headers
        headers = seq_ids.split('; ')
        
        species_set = set()
        for header in headers:
            species = parse_species_from_header(header)
            if species:
                species_set.add(species)
        
        exact_results.append({
            'Epitope_Sequence': motif,
            'Species_List': sorted(list(species_set)),
            'Total_Hits': row['UniqueSeqCount']
        })
    
    print()  # New line after progress
    
    exact_summary = pd.DataFrame(exact_results)
    exact_summary['Num_Species'] = exact_summary['Species_List'].apply(len)
    exact_summary['Species_Names'] = exact_summary['Species_List'].apply(
        lambda x: '; '.join(x) if len(x) > 0 else 'None'
    )
    exact_summary = exact_summary.sort_values('Num_Species', ascending=True)
    
    # Save detailed results
    exact_output = f"{output_dir}/exact_species_distribution.txt"
    exact_summary[['Epitope_Sequence', 'Num_Species', 'Total_Hits', 'Species_Names']].to_csv(
        exact_output, sep='\t', index=False
    )
    print(f"✓ Saved: {exact_output}")
    
    # Save detailed breakdown
    detailed_list = []
    for seq in exact_summary['Epitope_Sequence'].unique():
        seq_data = exact_summary[exact_summary['Epitope_Sequence'] == seq].iloc[0]
        
        for species in seq_data['Species_List']:
            detailed_list.append({
                'Epitope_Sequence': seq,
                'Species': species,
                'Total_Hits': seq_data['Total_Hits']
            })
    
    if len(detailed_list) > 0:
        detailed_df = pd.DataFrame(detailed_list)
        detailed_df = detailed_df.sort_values(['Epitope_Sequence', 'Species'])
        
        detailed_file = f"{output_dir}/exact_species_detailed.txt"
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        print(f"✓ Saved detailed breakdown: {detailed_file}")
    
    return exact_summary

def analyze_iedb_blast(output_dir):
    """Analyze IEDB BLAST with detailed output."""
    
    print("\n" + "="*70)
    print("ANALYZING IEDB BLAST (DETAILED)")
    print("="*70)
    
    iedb_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results/epitope_blast_COMBINED_fullheaders.txt'
    iedb_df = pd.read_csv(iedb_file, sep='\t')
    
    print(f"Total IEDB BLAST records: {len(iedb_df)}")
    
    # Parse species
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
    iedb_summary['Species_Names'] = iedb_summary['Species_List'].apply(
        lambda x: '; '.join(x) if len(x) > 0 else 'None'
    )
    iedb_summary = iedb_summary.sort_values('Num_Species', ascending=True)
    
    # Save
    iedb_output = f"{output_dir}/iedb_species_distribution.txt"
    iedb_summary[['Epitope_Sequence', 'Num_Species', 'Total_Hits', 'Species_Names']].to_csv(
        iedb_output, sep='\t', index=False
    )
    print(f"✓ Saved: {iedb_output}")
    
    # Save detailed breakdown
    iedb_anatomia_filtered = iedb_anatomia[iedb_anatomia['Species'].notna()].copy()
    
    detailed_list = []
    for seq in iedb_summary['Epitope_Sequence'].unique():
        seq_data = iedb_anatomia_filtered[iedb_anatomia_filtered['Epitope_Sequence'] == seq]
        
        for species in seq_data['Species'].unique():
            species_proteins = seq_data[seq_data['Species'] == species]
            
            detailed_list.append({
                'Epitope_Sequence': seq,
                'Species': species,
                'Num_Proteins': len(species_proteins),
                'Num_Unique_Proteins': species_proteins['Protein_Full_Header'].nunique(),
                'Example_Proteins': '; '.join(species_proteins['Protein_Full_Header'].unique()[:3].tolist())
            })
    
    if len(detailed_list) > 0:
        detailed_df = pd.DataFrame(detailed_list)
        detailed_df = detailed_df.sort_values(['Epitope_Sequence', 'Species'])
        
        detailed_file = f"{output_dir}/iedb_species_detailed.txt"
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        print(f"✓ Saved detailed breakdown: {detailed_file}")
    
    return iedb_summary

def main():
    """Main function to create all visualizations."""
    
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons'
    
    print("="*70)
    print("CREATING DETAILED VISUALIZATIONS FOR ALL METHODS")
    print("="*70)
    
    # Analyze and visualize Exact K-mer
    exact_df = analyze_exact_kmer(output_dir)
    visualize_single_method('Exact K-mer', exact_df, output_dir, color='#8da0cb')
    
    # Analyze and visualize IEDB BLAST
    iedb_df = analyze_iedb_blast(output_dir)
    visualize_single_method('IEDB BLAST', iedb_df, output_dir, color='#66c2a5')
    
    # Load and visualize Fuzzy (already analyzed)
    print("\n" + "="*70)
    print("RE-VISUALIZING FUZZY K-MER")
    print("="*70)
    
    fuzzy_file = f"{output_dir}/fuzzy_species_distribution.txt"
    fuzzy_df = pd.read_csv(fuzzy_file, sep='\t')
    visualize_single_method('Fuzzy K-mer', fuzzy_df, output_dir, color='#e78ac3')
    
    print("\n" + "="*70)
    print("ALL VISUALIZATIONS COMPLETE!")
    print("="*70)
    print(f"\nFiles created in: {output_dir}")
    print("  - exact_species_distribution.txt")
    print("  - exact_species_detailed.txt")
    print("  - exact_species_visualization.png")
    print("  - iedb_species_distribution.txt")
    print("  - iedb_species_detailed.txt")
    print("  - iedb_species_visualization.png")
    print("  - fuzzy_species_visualization.png")

if __name__ == "__main__":
    main()