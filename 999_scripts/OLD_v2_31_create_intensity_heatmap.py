import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import os
from datetime import datetime

def create_intensity_heatmap(metaproteome_dir, output_dir, log_dir):
    """
    Create publication-quality heatmap with intensity-based coloring and enhanced annotations.
    """
    
    print("="*80)
    print("CREATING INTENSITY-BASED VALIDATION HEATMAP")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load validated epitope data
    input_file = f"{metaproteome_dir}/REFINED_TIER1_with_metaproteome_validation.txt"
    print(f"\nLoading validated epitope data: {input_file}")
    
    df = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(df)} validated epitope clusters")
    
    # Load original intensity data for mapping
    intensity_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'
    intensity_df = pd.read_csv(intensity_file, sep='\t')
    
    # Create protein-to-intensity and localization mapping
    protein_intensity_map = dict(zip(intensity_df['Protein_ID'], intensity_df['Protein_Mean_Relative_Intensity']))
    protein_localization_map = dict(zip(intensity_df['Protein_ID'], intensity_df['Localization']))
    
    print(f"Loaded {len(protein_intensity_map)} protein intensity values")
    
    # Define species order (Gram+ then Gram-)
    species_order = [
        'B_longum', 'E_faecalis', 'L_casei', 'L_gasseri', 'L_reuteri', 
        'P_freudenreichii', 'S_salivarius', 'T_sanguinis',  # Gram+
        'A_muciniphila', 'B_fragilis', 'E_cloacae', 'F_nucleatum', 
        'F_prausnitzii', 'P_russellii', 'V_magna'  # Gram-
    ]
    
    # Process epitope data into heatmap format
    heatmap_data = []
    
    for _, epitope in df.iterrows():
        epitope_id = epitope['Cluster_ID']
        sequence = epitope['Representative_Sequence']
        proteins_str = epitope.get('Source_Proteins_All', '')
        
        # Initialize row data
        row_data = {
            'Epitope_ID': epitope_id,
            'Sequence': sequence,
            'Method': epitope['Representative_Method'],
            'Human_Identity': epitope['Max_Human_Identity'],
            'Total_Abundance': epitope.get('Total_Relative_Abundance', 0),
            'Metaproteome_Exact': epitope.get('Found_Exact_Match', False),
            'Metaproteome_BLAST': epitope.get('Found_BLAST_Match', False),
            'Surface_Accessible': False  # Will be determined from localization
        }
        
        # Initialize species intensities
        for species in species_order:
            row_data[species] = 0.0
        
        # Process source proteins
        if pd.notna(proteins_str) and proteins_str != '':
            surface_proteins = 0
            total_proteins = 0
            
            for protein_full in proteins_str.split(';'):
                protein_full = protein_full.strip()
                
                if protein_full.startswith('LITERATURE'):
                    # Handle literature data - extract species and mark as present (intensity = 1)
                    if 'L_casei' in protein_full.lower() or 'l_casei' in protein_full.lower():
                        row_data['L_casei'] = 1.0
                    elif 'L_reuteri' in protein_full.lower() or 'l_reuteri' in protein_full.lower():
                        row_data['L_reuteri'] = 1.0
                    elif 'E_faecalis' in protein_full.lower() or 'efaecalis' in protein_full.lower():
                        row_data['E_faecalis'] = 1.0
                    elif 'P_freudenreichii' in protein_full.lower():
                        row_data['P_freudenreichii'] = 1.0
                    elif 'B_longum' in protein_full.lower():
                        row_data['B_longum'] = 1.0
                    elif 'L_plantarum' in protein_full.lower():
                        row_data['L_plantarum'] = 1.0  # Note: not in species_order
                    
                elif protein_full.startswith('WP_'):
                    # Handle ANATOMIA data - extract protein ID and get intensity + localization
                    protein_id = protein_full.split()[0]
                    
                    if protein_id in protein_intensity_map:
                        intensity = protein_intensity_map[protein_id]
                        localization = protein_localization_map.get(protein_id, '')
                        
                        # Determine species from protein annotation
                        parts = protein_full.split('_')
                        if len(parts) >= 2:
                            species = f"{parts[0].split()[-1]}_{parts[1]}"
                            
                            if species in species_order and pd.notna(intensity):
                                # Add intensity to species (sum if multiple proteins)
                                row_data[species] += intensity
                        
                        # Check surface accessibility
                        total_proteins += 1
                        if localization in ['Cytoplasmic membrane', 'Outer membrane', 'Extracellular', 'Cellwall']:
                            surface_proteins += 1
            
            # Determine surface accessibility (majority of proteins)
            if total_proteins > 0:
                row_data['Surface_Accessible'] = (surface_proteins / total_proteins) > 0.5
        
        heatmap_data.append(row_data)
    
    # Convert to DataFrame
    heatmap_df = pd.DataFrame(heatmap_data)
    
    # Sort by total abundance (descending) then by species breadth
    species_counts = heatmap_df[species_order].apply(lambda x: (x > 0).sum(), axis=1)
    heatmap_df['Species_Count'] = species_counts
    heatmap_df = heatmap_df.sort_values(['Total_Abundance', 'Species_Count'], ascending=[False, False])
    
    # Take top 50 for visualization
    top_epitopes = heatmap_df.head(50).copy()
    
    print(f"\nCreating heatmap for top {len(top_epitopes)} epitopes")
    print(f"Species covered: {len(species_order)}")
    
    # Prepare data for heatmap
    intensity_matrix = top_epitopes[species_order].values
    epitope_labels = [f"{row['Sequence'][:25]}{'...' if len(row['Sequence']) > 25 else ''}" 
                     for _, row in top_epitopes.iterrows()]
    
    # Create figure
    fig, axes = plt.subplots(1, 1, figsize=(16, 20))
    
    # Create custom colormap for intensities
    # For ANATOMIA data: white (0) -> light blue (low) -> dark blue (high)
    # For Literature data: distinct color (orange)
    
    # Normalize intensity values for better visualization
    intensity_matrix_norm = intensity_matrix.copy()
    
    # Separate literature (1.0) from ANATOMIA intensities
    anatomia_mask = (intensity_matrix != 1.0) & (intensity_matrix > 0)
    literature_mask = intensity_matrix == 1.0
    
    # Normalize ANATOMIA intensities to 0.1-0.9 range (keeping 1.0 for literature)
    if np.any(anatomia_mask):
        anatomia_values = intensity_matrix[anatomia_mask]
        if len(anatomia_values) > 0:
            # Log transform for better visualization
            log_values = np.log10(anatomia_values + 1e-6)
            min_log = log_values.min()
            max_log = log_values.max()
            
            if max_log > min_log:
                normalized_log = 0.1 + 0.8 * (log_values - min_log) / (max_log - min_log)
                intensity_matrix_norm[anatomia_mask] = normalized_log
    
    # Create custom colormap
    colors = ['white', '#e6f3ff', '#b3d9ff', '#4da6ff', '#0066cc', '#003d7a']  # White to dark blue
    literature_color = '#ff7f0e'  # Orange for literature
    
    cmap = LinearSegmentedColormap.from_list('intensity', colors, N=100)
    
    # Create the heatmap
    im = axes.imshow(intensity_matrix_norm, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    
    # Overlay literature data with distinct color
    literature_coords = np.where(literature_mask)
    for i, j in zip(literature_coords[0], literature_coords[1]):
        axes.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, 
                                   facecolor=literature_color, alpha=0.8, edgecolor='none'))
    
    # Set labels
    axes.set_xticks(range(len(species_order)))
    axes.set_xticklabels([s.replace('_', ' ') for s in species_order], rotation=45, ha='right')
    axes.set_yticks(range(len(epitope_labels)))
    axes.set_yticklabels(epitope_labels, fontsize=8)
    
    # Add species group separation
    gram_pos_count = 8  # Number of Gram+ species
    axes.axvline(x=gram_pos_count-0.5, color='black', linewidth=2, alpha=0.7)
    
    # Add annotation columns
    # Calculate positions for annotation columns
    n_species = len(species_order)
    annotation_start = n_species + 0.5
    
    # Metaproteome detection (split into exact and BLAST)
    meta_exact_col = annotation_start
    meta_blast_col = annotation_start + 0.5
    
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        # Exact match
        if row['Metaproteome_Exact']:
            axes.add_patch(plt.Rectangle((meta_exact_col-0.25, i-0.5), 0.5, 1, 
                                       facecolor='#2ca02c', alpha=0.8))  # Green
        
        # BLAST match (only if no exact match)
        if row['Metaproteome_BLAST'] and not row['Metaproteome_Exact']:
            axes.add_patch(plt.Rectangle((meta_blast_col-0.25, i-0.5), 0.5, 1, 
                                       facecolor='#90EE90', alpha=0.8))  # Light green
    
    # Method annotation
    method_col = annotation_start + 1.5
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        method = row['Method']
        if 'IEDB' in method:
            color = '#9467bd'  # Purple
        else:
            color = '#ff7f0e'  # Orange
        
        axes.add_patch(plt.Rectangle((method_col-0.5, i-0.5), 1, 1, 
                                   facecolor=color, alpha=0.8))
    
    # Surface accessibility (Fish)
    fish_col = annotation_start + 2.5
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        if row['Surface_Accessible']:
            axes.add_patch(plt.Rectangle((fish_col-0.5, i-0.5), 1, 1, 
                                       facecolor='#d62728', alpha=0.8))  # Red
        else:
            axes.add_patch(plt.Rectangle((fish_col-0.5, i-0.5), 1, 1, 
                                       facecolor='#7f7f7f', alpha=0.8))  # Gray
    
    # Set axis limits
    axes.set_xlim(-0.5, fish_col + 0.5)
    axes.set_ylim(len(epitope_labels)-0.5, -0.5)
    
    # Add column headers for annotations
    axes.text(meta_exact_col, -1, 'Meta\nExact', ha='center', va='bottom', fontsize=10, fontweight='bold')
    axes.text(meta_blast_col, -1, 'Meta\nBLAST', ha='center', va='bottom', fontsize=10, fontweight='bold')
    axes.text(method_col, -1, 'Method', ha='center', va='bottom', fontsize=10, fontweight='bold')
    axes.text(fish_col, -1, 'Fish', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Add legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', label='Not detected'),
        plt.Rectangle((0,0),1,1, facecolor='#4da6ff', label='ANATOMIA intensity'),
        plt.Rectangle((0,0),1,1, facecolor=literature_color, label='Literature data'),
        plt.Rectangle((0,0),1,1, facecolor='#2ca02c', label='Metaproteome: Exact'),
        plt.Rectangle((0,0),1,1, facecolor='#90EE90', label='Metaproteome: BLAST'),
        plt.Rectangle((0,0),1,1, facecolor='#9467bd', label='Method: IEDB'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', label='Method: K-mer'),
        plt.Rectangle((0,0),1,1, facecolor='#d62728', label='Surface accessible'),
        plt.Rectangle((0,0),1,1, facecolor='#7f7f7f', label='Not surface accessible'),
    ]
    
    axes.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=10)
    
    # Set title
    plt.title(f'TOP{len(top_epitopes)} Epitopes: Multi-Source Validation with ANATOMIA Intensities\n' +
             f'Species Coverage (Gram+ | Gram-)', 
             fontsize=16, fontweight='bold', pad=20)
    
    # Add Gram+/Gram- labels
    axes.text(gram_pos_count/2 - 0.5, -2.5, 'Gram+', ha='center', va='center', 
             fontsize=12, fontweight='bold')
    axes.text(gram_pos_count + (len(species_order)-gram_pos_count)/2 - 0.5, -2.5, 'Gram-', 
             ha='center', va='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    output_file = f"{output_dir}/intensity_heatmap_multisource_validation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved heatmap: {output_file}")
    plt.close()
    
    # Save data used for heatmap
    data_file = f"{output_dir}/heatmap_data_top{len(top_epitopes)}.txt"
    top_epitopes.to_csv(data_file, sep='\t', index=False)
    print(f"✓ Saved data: {data_file}")
    
    # Print summary
    print(f"\n" + "="*80)
    print("HEATMAP SUMMARY")
    print("="*80)
    print(f"Epitopes shown: {len(top_epitopes)}")
    print(f"Species covered: {len(species_order)}")
    print(f"Metaproteome exact matches: {top_epitopes['Metaproteome_Exact'].sum()}")
    print(f"Metaproteome BLAST matches: {top_epitopes['Metaproteome_BLAST'].sum()}")
    print(f"Surface accessible: {top_epitopes['Surface_Accessible'].sum()}")
    print(f"IEDB-derived: {top_epitopes['Method'].str.contains('IEDB').sum()}")
    print(f"Mean abundance: {top_epitopes['Total_Abundance'].mean():.4f}")

if __name__ == "__main__":
    metaproteome_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation'
    output_dir = f"{metaproteome_dir}/visualizations"
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    create_intensity_heatmap(metaproteome_dir, output_dir, log_dir)