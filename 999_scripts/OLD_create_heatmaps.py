import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Circle
import os

def create_epitope_species_heatmap(epitope_file, intensity_file, metaproteome_file, output_dir):
    """
    Create heatmap of epitope presence across species with annotations.
    Colorblind-friendly using viridis colormap.
    """
    
    print("="*70)
    print("CREATING EPITOPE × SPECIES HEATMAP")
    print("="*70)
    
    # Load data
    print("\nLoading data...")
    epitope_df = pd.read_csv(epitope_file, sep='\t')
    intensity_df = pd.read_csv(intensity_file, sep='\t')
    metaproteome_df = pd.read_csv(metaproteome_file, sep='\t')
    
    print(f"Epitopes: {len(epitope_df)}")
    
    # Load protein mappings to get epitope-species-protein relationships
    mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/TOP50_protein_mapping.txt'
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Parse species from Full_Protein_Header for ANATOMIA proteins
    def parse_species_from_header(header, data_source):
        if pd.isna(header) or data_source != 'ANATOMIA':
            return None
        
        header_str = str(header)
        if ' ' not in header_str:
            return None
        
        # Format: ProteinID Species_OG_Conservation_Localization
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
    
    mapping_df['Parsed_Species'] = mapping_df.apply(
        lambda row: parse_species_from_header(row['Full_Protein_Header'], row['Data_Source']),
        axis=1
    )
    
    # Filter to ANATOMIA only
    anatomia_mapping = mapping_df[mapping_df['Data_Source'] == 'ANATOMIA'].copy()
    
    # Get unique species
    species_list = sorted(anatomia_mapping['Parsed_Species'].dropna().unique())
    print(f"Species: {len(species_list)}")
    
    # Create epitope × species matrix for hit counts
    epitope_ids = epitope_df['Cluster_ID'].tolist()
    hit_matrix = pd.DataFrame(0, index=epitope_ids, columns=species_list)
    
    # Fill matrix with hit counts
    for epitope_id in epitope_ids:
        epitope_proteins = anatomia_mapping[anatomia_mapping['Epitope_ID'] == epitope_id]
        species_counts = epitope_proteins['Parsed_Species'].value_counts()
        
        for species, count in species_counts.items():
            if species in hit_matrix.columns:
                hit_matrix.loc[epitope_id, species] = count
    
    # Create metaproteome presence matrix (binary: has intensity > 0)
    metaproteome_matrix = pd.DataFrame(0, index=epitope_ids, columns=species_list)
    
    for _, row in metaproteome_df.iterrows():
        epitope_id = row['Cluster_ID']
        if epitope_id in metaproteome_matrix.index and row['Num_Metaproteome_Proteins'] > 0:
            # Mark all species where this epitope was found
            # Note: metaproteome doesn't have species breakdown, so this is approximate
            metaproteome_matrix.loc[epitope_id, :] = 1  # Mark presence
    
    # Sort epitopes by Priority Score (descending)
    epitope_order = epitope_df.sort_values('Priority_Score', ascending=False)['Cluster_ID'].tolist()
    hit_matrix = hit_matrix.loc[epitope_order]
    metaproteome_matrix = metaproteome_matrix.loc[epitope_order]
    
    # Prepare annotations
    epitope_annotations = epitope_df.set_index('Cluster_ID').loc[epitope_order]
    
    # Create figure
    fig = plt.figure(figsize=(16, 12))
    
    # Define grid
    gs = fig.add_gridspec(1, 4, width_ratios=[0.3, 0.1, 3, 0.2], wspace=0.05)
    
    # Left annotations
    ax_method = fig.add_subplot(gs[0, 0])
    ax_fish = fig.add_subplot(gs[0, 1])
    ax_heatmap = fig.add_subplot(gs[0, 2])
    ax_colorbar = fig.add_subplot(gs[0, 3])
    
    # Method annotation
    method_colors = {
        'HIGH_CONFIDENCE': '#E69F00',  # Orange (colorblind-safe)
        'IEDB': '#56B4E9',              # Sky blue
        'Kmer_EXACT': '#009E73',        # Green
        'Kmer_FUZZY': '#CC79A7'         # Pink
    }
    
    def get_method_color(methods_str):
        if 'HIGH_CONFIDENCE' in str(methods_str):
            return method_colors['HIGH_CONFIDENCE']
        elif 'IEDB' in str(methods_str):
            return method_colors['IEDB']
        elif 'Kmer_EXACT' in str(methods_str):
            return method_colors['Kmer_EXACT']
        else:
            return method_colors['Kmer_FUZZY']
    
    method_annotation = epitope_annotations['Methods_Used'].apply(get_method_color)
    
    for i, color in enumerate(method_annotation):
        ax_method.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5))
    
    ax_method.set_xlim(0, 1)
    ax_method.set_ylim(0, len(epitope_order))
    ax_method.set_yticks([])
    ax_method.set_xticks([])
    ax_method.set_ylabel('')
    ax_method.set_title('Method', fontsize=10)
    ax_method.invert_yaxis()
    
    # Fishability annotation
    fish_colors = {'Yes': '#0072B2', 'No': '#D55E00'}  # Blue/Orange colorblind-safe
    fishability_annotation = epitope_annotations['Suitable_For_Fishing'].map(fish_colors)
    
    for i, color in enumerate(fishability_annotation):
        ax_fish.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='white', linewidth=0.5))
    
    ax_fish.set_xlim(0, 1)
    ax_fish.set_ylim(0, len(epitope_order))
    ax_fish.set_yticks([])
    ax_fish.set_xticks([])
    ax_fish.set_title('Fish', fontsize=10, rotation=0)
    ax_fish.invert_yaxis()
    
    # Main heatmap
    sns.heatmap(
        hit_matrix,
        ax=ax_heatmap,
        cmap='viridis',  # Colorblind-friendly
        cbar_ax=ax_colorbar,
        linewidths=0.5,
        linecolor='lightgray',
        square=False,
        cbar_kws={'label': 'Number of ANATOMIA hits'}
    )
    
    # Add dots for metaproteome presence
    # (This is complex - for now skip, or add in next version)
    
    # Labels
    ax_heatmap.set_xlabel('Species', fontsize=12, fontweight='bold')
    ax_heatmap.set_ylabel('Epitope (ordered by Priority Score)', fontsize=12, fontweight='bold')
    ax_heatmap.set_yticklabels(epitope_annotations['Representative_Sequence'], fontsize=8)
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    
    plt.suptitle('TOP50 Epitopes: Presence Across Species', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    # Add legend for annotations
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=method_colors['HIGH_CONFIDENCE'], label='HIGH_CONFIDENCE'),
        Patch(facecolor=method_colors['IEDB'], label='IEDB'),
        Patch(facecolor=method_colors['Kmer_EXACT'], label='Kmer_EXACT'),
        Patch(facecolor=method_colors['Kmer_FUZZY'], label='Kmer_FUZZY'),
        Patch(facecolor=fish_colors['Yes'], label='Fishable'),
        Patch(facecolor=fish_colors['No'], label='Not Fishable')
    ]
    
    fig.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.95), 
               fontsize=9, frameon=True, title='Annotations')
    
    # Save
    output_file = f"{output_dir}/heatmap_epitope_species_TOP50.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    
    plt.close()
    
    # Create summary statistics
    print("\nHeatmap statistics:")
    print(f"  Total cells: {hit_matrix.size}")
    print(f"  Cells with hits: {(hit_matrix > 0).sum().sum()}")
    print(f"  Max hits in single cell: {hit_matrix.max().max()}")
    print(f"  Species with most epitopes: {(hit_matrix > 0).sum().idxmax()} ({(hit_matrix > 0).sum().max()} epitopes)")
    
    return hit_matrix

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/007_visualizations'
os.makedirs(output_dir, exist_ok=True)

epitope_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt'
intensity_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_intensity_analysis/TOP50_FINAL_with_intensities.txt'
metaproteome_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation/TOP50_FINAL_metaproteome_search.txt'

hit_matrix = create_epitope_species_heatmap(epitope_file, intensity_file, metaproteome_file, output_dir)

print("\n" + "="*70)
print("HEATMAP CREATION COMPLETE")
print("="*70)