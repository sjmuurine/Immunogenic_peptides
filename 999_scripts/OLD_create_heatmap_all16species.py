import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon
import os

def create_multi_source_heatmap(output_dir):
    """
    Heatmap showing ALL 16 species, even if some have no TOP50 epitopes.
    """
    
    print("="*70)
    print("CREATING MULTI-SOURCE EPITOPE × SPECIES HEATMAP (ALL 16 SPECIES)")
    print("="*70)
    
    # Load data
    print("\nLoading data...")
    epitope_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt'
    epitope_df = pd.read_csv(epitope_file, sep='\t')
    
    detailed_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_cluster_detailed_v2.txt'
    detailed_df = pd.read_csv(detailed_file, sep='\t')
    
    cluster_to_epitopes = {}
    for cluster_id in detailed_df['Cluster_ID'].unique():
        epitope_ids = detailed_df[detailed_df['Cluster_ID'] == cluster_id]['Epitope_ID'].tolist()
        cluster_to_epitopes[cluster_id] = epitope_ids
    
    mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/TOP50_protein_mapping.txt'
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    def parse_species_from_header(header, data_source):
        if pd.isna(header):
            return None
        header_str = str(header)
        if data_source == 'LITERATURE':
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
            return None
        if ' ' not in header_str:
            return None
        parts = header_str.split(' ', 1)[1].split('_')
        og_idx = next((i for i, p in enumerate(parts) if p.startswith('OG') or p == 'unique'), None)
        return '_'.join(parts[:og_idx]) if og_idx and og_idx > 0 else None
    
    mapping_df['Parsed_Species'] = mapping_df.apply(
        lambda row: parse_species_from_header(row['Full_Protein_Header'], row['Data_Source']), axis=1
    )
    
    # Define ALL 16 species in Gram order
    gram_positive = ['B_longum', 'E_faecalis', 'L_casei', 'L_plantarum', 'L_reuteri', 
                     'P_russellii', 'S_salivarius', 'T_sanguinis']
    gram_negative = ['A_muciniphila', 'B_fragilis', 'E_cloacae', 'F_nucleatum', 
                     'F_prausnitzii', 'L_fermentum', 'L_gasseri', 'V_magna']
    
    # Complete species list (all 16)
    species_list = gram_positive + gram_negative
    
    print(f"\nAll 16 species (Gram+ then Gram-):")
    print(f"  Gram+: {gram_positive}")
    print(f"  Gram-: {gram_negative}")
    
    epitope_order = epitope_df.sort_values('Priority_Score', ascending=False)['Cluster_ID'].tolist()
    
    # Create matrices for all 16 species
    anatomia_matrix = pd.DataFrame(0, index=epitope_order, columns=species_list)
    literature_matrix = pd.DataFrame(0, index=epitope_order, columns=species_list)
    
    for cluster_id in epitope_order:
        for epitope_id in cluster_to_epitopes.get(cluster_id, []):
            epitope_matches = mapping_df[mapping_df['Epitope_ID'] == epitope_id]
            for species in epitope_matches[epitope_matches['Data_Source'] == 'ANATOMIA']['Parsed_Species'].dropna():
                if species in anatomia_matrix.columns:
                    anatomia_matrix.loc[cluster_id, species] = 1
            for species in epitope_matches[epitope_matches['Data_Source'] == 'LITERATURE']['Parsed_Species'].dropna():
                if species in literature_matrix.columns:
                    literature_matrix.loc[cluster_id, species] = 1
    
    metaproteome_df = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation/TOP50_FINAL_metaproteome_search.txt', sep='\t')
    metaproteome_present = set(metaproteome_df[metaproteome_df['Found_In_Metaproteome'] == 'Yes']['Cluster_ID'])
    
    epitope_annotations = epitope_df.set_index('Cluster_ID').loc[epitope_order]
    
    # Colors
    method_colors = {'HIGH_CONFIDENCE': '#E69F00', 'IEDB': '#56B4E9', 'Kmer_EXACT': '#009E73', 'Kmer_FUZZY': '#CC79A7'}
    fish_colors = {'Yes': '#0072B2', 'No': '#D55E00'}
    metaproteome_color = '#9b59b6'
    
    def get_method_color(m):
        if pd.isna(m): return '#999999'
        m = str(m)
        return next((method_colors[k] for k in ['HIGH_CONFIDENCE', 'IEDB', 'Kmer_EXACT'] if k in m), method_colors['Kmer_FUZZY'])
    
    # Create figure
    fig = plt.figure(figsize=(24, 14))  # Wider for 16 species
    
    ax_sequences = plt.axes([-0.08, 0.08, 0.23, 0.85])
    ax_heatmap = plt.axes([0.21, 0.08, 0.53, 0.85])  # Slightly wider
    ax_metaproteome = plt.axes([0.75, 0.08, 0.02, 0.85])
    ax_method = plt.axes([0.78, 0.08, 0.04, 0.85])
    ax_fish = plt.axes([0.83, 0.08, 0.02, 0.85])
    
    # Sequences
    ax_sequences.set_xlim(0, 1)
    ax_sequences.set_ylim(0, len(epitope_order))
    ax_sequences.invert_yaxis()
    ax_sequences.set_yticks(np.arange(len(epitope_order)) + 0.5)
    ax_sequences.set_yticklabels(epitope_annotations['Representative_Sequence'], fontsize=7, family='monospace')
    ax_sequences.yaxis.tick_right()
    ax_sequences.set_xticks([])
    ax_sequences.tick_params(left=False, right=False)
    for spine in ax_sequences.spines.values():
        spine.set_visible(False)
    
    # Heatmap
    ax_heatmap.set_xlim(0, len(species_list))
    ax_heatmap.set_ylim(0, len(epitope_order))
    ax_heatmap.invert_yaxis()
    
    for i, epitope_id in enumerate(epitope_order):
        for j, species in enumerate(species_list):
            ax_heatmap.add_patch(plt.Rectangle((j, i), 1, 1, facecolor='white', edgecolor='lightgray', linewidth=0.5))
            
            has_anat = anatomia_matrix.loc[epitope_id, species] > 0
            has_lit = literature_matrix.loc[epitope_id, species] > 0
            
            if has_anat:
                ax_heatmap.add_patch(plt.Rectangle((j+0.15, i+0.15), 0.7, 0.7, facecolor='#2a9d8f', edgecolor='none'))
            
            if has_lit:
                ax_heatmap.add_patch(Polygon([[j+0.7, i+0.1], [j+0.9, i+0.1], [j+0.9, i+0.3]], facecolor='#e76f51', edgecolor='none'))
    
    # Add visual separator between Gram+ and Gram-
    gram_pos_count = len(gram_positive)
    ax_heatmap.axvline(x=gram_pos_count, color='black', linewidth=2, linestyle='--', alpha=0.5)
    
    ax_heatmap.set_xticks(np.arange(len(species_list)) + 0.5)
    ax_heatmap.set_xticklabels(species_list, rotation=45, ha='right', fontsize=9)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel('Species (Gram+ | Gram-)', fontsize=12, fontweight='bold')
    ax_heatmap.tick_params(left=False, bottom=False)
    
    # Metaproteome
    ax_metaproteome.set_xlim(0, 1)
    ax_metaproteome.set_ylim(0, len(epitope_order))
    ax_metaproteome.invert_yaxis()
    
    for i, epitope_id in enumerate(epitope_order):
        color = metaproteome_color if epitope_id in metaproteome_present else 'white'
        ax_metaproteome.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='lightgray', linewidth=0.5))
    
    ax_metaproteome.set_yticks([])
    ax_metaproteome.set_xticks([])
    ax_metaproteome.set_title('Meta', fontsize=10, fontweight='bold', pad=5)
    for spine in ax_metaproteome.spines.values():
        spine.set_visible(False)
    
    # Method
    ax_method.set_xlim(0, 1)
    ax_method.set_ylim(0, len(epitope_order))
    ax_method.invert_yaxis()
    
    for i, methods in enumerate(epitope_annotations['Methods_Used']):
        ax_method.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=get_method_color(methods), edgecolor='white', linewidth=0.5))
    
    ax_method.set_yticks([])
    ax_method.set_xticks([])
    ax_method.set_title('Method', fontsize=10, fontweight='bold', pad=5)
    for spine in ax_method.spines.values():
        spine.set_visible(False)
    
    # Fish
    ax_fish.set_xlim(0, 1)
    ax_fish.set_ylim(0, len(epitope_order))
    ax_fish.invert_yaxis()
    
    for i, fishable in enumerate(epitope_annotations['Suitable_For_Fishing']):
        ax_fish.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=fish_colors.get(fishable, '#CCC'), edgecolor='white', linewidth=0.5))
    
    ax_fish.set_yticks([])
    ax_fish.set_xticks([])
    ax_fish.set_title('Fish', fontsize=10, fontweight='bold', pad=5)
    for spine in ax_fish.spines.values():
        spine.set_visible(False)
    
    # Title
    fig.text(0.5, 0.97, 'TOP50 Epitopes: Multi-Source Validation Across All 16 Species', 
             ha='center', fontsize=16, fontweight='bold')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=method_colors['HIGH_CONFIDENCE'], label='HIGH_CONFIDENCE', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['IEDB'], label='IEDB', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['Kmer_EXACT'], label='Kmer_EXACT', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['Kmer_FUZZY'], label='Kmer_FUZZY', edgecolor='black'),
        mpatches.Patch(facecolor='white', label='', edgecolor='none'),
        mpatches.Patch(facecolor=fish_colors['Yes'], label='Fishable', edgecolor='black'),
        mpatches.Patch(facecolor=fish_colors['No'], label='Not Fishable', edgecolor='black'),
        mpatches.Patch(facecolor='white', label='', edgecolor='none'),
        mpatches.Patch(facecolor='#2a9d8f', label='ANATOMIA ■', edgecolor='black'),
        mpatches.Patch(facecolor='#e76f51', label='LITERATURE ▲', edgecolor='black'),
        mpatches.Patch(facecolor=metaproteome_color, label='Metaproteome', edgecolor='black'),
    ]
    
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.92, 0.96), 
               fontsize=9, frameon=True, ncol=1, title='Legend', title_fontsize=10)
    
    output_file = f"{output_dir}/heatmap_multisource_TOP50_ALL16species.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Saved: {output_file}")
    plt.close()
    
    # Print species coverage summary
    print("\n" + "="*70)
    print("SPECIES COVERAGE IN TOP50")
    print("="*70)
    
    for species in species_list:
        count = (anatomia_matrix[species] + literature_matrix[species] > 0).sum()
        print(f"{species:20} {count:2} epitopes")
    
    print("\nCOMPLETE!")

output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/007_visualizations'
os.makedirs(output_dir, exist_ok=True)
create_multi_source_heatmap(output_dir)