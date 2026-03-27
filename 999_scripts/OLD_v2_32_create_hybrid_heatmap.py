import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import os
from datetime import datetime

def create_hybrid_heatmap(metaproteome_dir, output_dir, log_dir):
    """
    Create hybrid heatmap with binary species presence and abundance as separate annotation.
    """
    
    print("="*80)
    print("CREATING HYBRID HEATMAP (BINARY + ABUNDANCE)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load validated epitope data
    input_file = f"{metaproteome_dir}/REFINED_TIER1_with_metaproteome_validation.txt"
    print(f"\nLoading validated epitope data: {input_file}")
    
    df = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(df)} validated epitope clusters")
    
    # Load original intensity data for localization mapping
    intensity_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'
    intensity_df = pd.read_csv(intensity_file, sep='\t')
    protein_localization_map = dict(zip(intensity_df['Protein_ID'], intensity_df['Localization']))
    
    print(f"Loaded localization data for {len(protein_localization_map)} proteins")
    
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
            'Surface_Accessible': False
        }
        
        # Initialize species presence (binary)
        for species in species_order:
            row_data[species] = 0  # Binary: 0 = absent, 1 = present
        
        # Process source proteins for species presence and surface accessibility
        if pd.notna(proteins_str) and proteins_str != '':
            surface_proteins = 0
            total_proteins = 0
            
            for protein_full in proteins_str.split(';'):
                protein_full = protein_full.strip()
                
                if protein_full.startswith('LITERATURE'):
                    # Handle literature data - mark species as present
                    if 'L_casei' in protein_full.lower() or 'l_casei' in protein_full.lower():
                        row_data['L_casei'] = 1
                    elif 'L_reuteri' in protein_full.lower() or 'l_reuteri' in protein_full.lower():
                        row_data['L_reuteri'] = 1
                    elif 'E_faecalis' in protein_full.lower() or 'efaecalis' in protein_full.lower():
                        row_data['E_faecalis'] = 1
                    elif 'P_freudenreichii' in protein_full.lower():
                        row_data['P_freudenreichii'] = 1
                    elif 'B_longum' in protein_full.lower():
                        row_data['B_longum'] = 1
                    elif 'L_plantarum' in protein_full.lower():
                        pass  # L_plantarum not in our species_order
                    
                elif protein_full.startswith('WP_'):
                    # Handle ANATOMIA data
                    protein_id = protein_full.split()[0]
                    localization = protein_localization_map.get(protein_id, '')
                    
                    # Determine species from protein annotation
                    parts = protein_full.split('_')
                    if len(parts) >= 2:
                        species = f"{parts[0].split()[-1]}_{parts[1]}"
                        if species in species_order:
                            row_data[species] = 1  # Mark as present
                    
                    # Check surface accessibility
                    total_proteins += 1
                    if localization in ['Cytoplasmic membrane', 'Outer membrane', 'Extracellular', 'Cellwall']:
                        surface_proteins += 1
            
            # Determine overall surface accessibility
            if total_proteins > 0:
                row_data['Surface_Accessible'] = (surface_proteins / total_proteins) > 0.5
        
        heatmap_data.append(row_data)
    
    # Convert to DataFrame
    heatmap_df = pd.DataFrame(heatmap_data)
    
    # Sort by total abundance (descending) then by species breadth
    species_counts = heatmap_df[species_order].sum(axis=1)
    heatmap_df['Species_Count'] = species_counts
    heatmap_df = heatmap_df.sort_values(['Total_Abundance', 'Species_Count'], ascending=[False, False])
    
    # Take top 50 for visualization
    top_epitopes = heatmap_df.head(50).copy()
    
    print(f"\nCreating hybrid heatmap for top {len(top_epitopes)} epitopes")
    print(f"Species covered: {len(species_order)}")
    print(f"Abundance range: {top_epitopes['Total_Abundance'].min():.4f} - {top_epitopes['Total_Abundance'].max():.4f}")
    
    # Prepare data for heatmap
    species_matrix = top_epitopes[species_order].values
    epitope_labels = [f"{row['Sequence'][:25]}{'...' if len(row['Sequence']) > 25 else ''}" 
                     for _, row in top_epitopes.iterrows()]
    
    # Create figure
    fig = plt.figure(figsize=(18, 22))
    
    # Define layout: species heatmap + annotation columns
    # Species columns + abundance column + metaproteome columns + method + fish
    n_species = len(species_order)
    total_width = n_species + 5  # 5 annotation columns
    
    # Create main heatmap
    main_ax = plt.subplot2grid((1, total_width), (0, 0), colspan=n_species, fig=fig)
    
    # Plot binary species presence
    im = main_ax.imshow(species_matrix, cmap='RdBu_r', aspect='auto', vmin=0, vmax=1)
    
    # Customize species heatmap
    main_ax.set_xticks(range(len(species_order)))
    main_ax.set_xticklabels([s.replace('_', ' ') for s in species_order], rotation=45, ha='right')
    main_ax.set_yticks(range(len(epitope_labels)))
    main_ax.set_yticklabels(epitope_labels, fontsize=9)
    
    # Add Gram+/Gram- separation
    gram_pos_count = 8
    main_ax.axvline(x=gram_pos_count-0.5, color='black', linewidth=3, alpha=0.8)
    
    # Set limits
    main_ax.set_xlim(-0.5, n_species-0.5)
    main_ax.set_ylim(len(epitope_labels)-0.5, -0.5)
    
    # Add Gram+/Gram- labels
    main_ax.text(gram_pos_count/2 - 0.5, -2, 'Gram+', ha='center', va='center', 
                fontsize=12, fontweight='bold')
    main_ax.text(gram_pos_count + (n_species-gram_pos_count)/2 - 0.5, -2, 'Gram-', 
                ha='center', va='center', fontsize=12, fontweight='bold')
    
    # Create annotation columns
    col_positions = [n_species + i + 0.5 for i in range(5)]
    col_names = ['Abundance', 'Meta\nExact', 'Meta\nBLAST', 'Method', 'Fish']
    
    # 1. Abundance column (color gradient)
    abundance_ax = plt.subplot2grid((1, total_width), (0, n_species), fig=fig)
    
    # Create abundance colormap
    abundance_values = top_epitopes['Total_Abundance'].values
    abundance_normalized = (abundance_values - abundance_values.min()) / (abundance_values.max() - abundance_values.min() + 1e-10)
    
    abundance_cmap = LinearSegmentedColormap.from_list('abundance', ['white', '#08519c'], N=100)
    abundance_colors = abundance_cmap(abundance_normalized)
    
    # Plot abundance as colored rectangles
    for i, color in enumerate(abundance_colors):
        abundance_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='lightgray', linewidth=0.5))
    
    abundance_ax.set_xlim(0, 1)
    abundance_ax.set_ylim(0, len(top_epitopes))
    abundance_ax.set_xticks([])
    abundance_ax.set_yticks([])
    abundance_ax.set_title('Abundance', fontsize=10, fontweight='bold', pad=10)
    
    # 2. Metaproteome Exact column
    meta_exact_ax = plt.subplot2grid((1, total_width), (0, n_species + 1), fig=fig)
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        color = '#2ca02c' if row['Metaproteome_Exact'] else 'white'
        meta_exact_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, 
                                            edgecolor='lightgray', linewidth=0.5))
    
    meta_exact_ax.set_xlim(0, 1)
    meta_exact_ax.set_ylim(0, len(top_epitopes))
    meta_exact_ax.set_xticks([])
    meta_exact_ax.set_yticks([])
    meta_exact_ax.set_title('Meta\nExact', fontsize=10, fontweight='bold', pad=10)
    
    # 3. Metaproteome BLAST column
    meta_blast_ax = plt.subplot2grid((1, total_width), (0, n_species + 2), fig=fig)
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        # Only show BLAST if no exact match
        if row['Metaproteome_BLAST'] and not row['Metaproteome_Exact']:
            color = '#90EE90'  # Light green
        else:
            color = 'white'
        meta_blast_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color,
                                            edgecolor='lightgray', linewidth=0.5))
    
    meta_blast_ax.set_xlim(0, 1)
    meta_blast_ax.set_ylim(0, len(top_epitopes))
    meta_blast_ax.set_xticks([])
    meta_blast_ax.set_yticks([])
    meta_blast_ax.set_title('Meta\nBLAST', fontsize=10, fontweight='bold', pad=10)
    
    # 4. Method column
    method_ax = plt.subplot2grid((1, total_width), (0, n_species + 3), fig=fig)
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        method = row['Method']
        color = '#9467bd' if 'IEDB' in method else '#ff7f0e'  # Purple for IEDB, Orange for K-mer
        method_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color,
                                        edgecolor='lightgray', linewidth=0.5))
    
    method_ax.set_xlim(0, 1)
    method_ax.set_ylim(0, len(top_epitopes))
    method_ax.set_xticks([])
    method_ax.set_yticks([])
    method_ax.set_title('Method', fontsize=10, fontweight='bold', pad=10)
    
    # 5. Surface accessibility (Fish) column
    fish_ax = plt.subplot2grid((1, total_width), (0, n_species + 4), fig=fig)
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        color = '#d62728' if row['Surface_Accessible'] else '#7f7f7f'  # Red vs Gray
        fish_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color,
                                      edgecolor='lightgray', linewidth=0.5))
    
    fish_ax.set_xlim(0, 1)
    fish_ax.set_ylim(0, len(top_epitopes))
    fish_ax.set_xticks([])
    fish_ax.set_yticks([])
    fish_ax.set_title('Fish', fontsize=10, fontweight='bold', pad=10)
    
    # Create legend
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', label='Not present'),
        plt.Rectangle((0,0),1,1, facecolor='#b2182b', label='Present'),
        plt.Rectangle((0,0),1,1, facecolor='#08519c', label='High abundance'),
        plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='gray', label='Low abundance'),
        plt.Rectangle((0,0),1,1, facecolor='#2ca02c', label='Metaproteome: Exact'),
        plt.Rectangle((0,0),1,1, facecolor='#90EE90', label='Metaproteome: BLAST only'),
        plt.Rectangle((0,0),1,1, facecolor='#9467bd', label='Method: IEDB'),
        plt.Rectangle((0,0),1,1, facecolor='#ff7f0e', label='Method: K-mer'),
        plt.Rectangle((0,0),1,1, facecolor='#d62728', label='Surface accessible'),
        plt.Rectangle((0,0),1,1, facecolor='#7f7f7f', label='Not surface accessible'),
    ]
    
    # Position legend to the right of all columns
    fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.92, 0.5), fontsize=11)
    
    # Set main title
    fig.suptitle(f'TOP{len(top_epitopes)} Epitopes: Multi-Source Validation with Binary Presence + Abundance\n' +
                f'Species Coverage (Gram+ | Gram-)', 
                fontsize=16, fontweight='bold', y=0.95)
    
    # Adjust layout
    plt.subplots_adjust(left=0.3, right=0.85, top=0.9, bottom=0.1, wspace=0.1)
    
    # Save figure
    output_file = f"{output_dir}/hybrid_heatmap_binary_abundance.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved hybrid heatmap: {output_file}")
    plt.close()
    
    # Save data
    data_file = f"{output_dir}/hybrid_heatmap_data_top{len(top_epitopes)}.txt"
    top_epitopes.to_csv(data_file, sep='\t', index=False)
    print(f"✓ Saved data: {data_file}")
    
    # Print summary
    print(f"\n" + "="*80)
    print("HYBRID HEATMAP SUMMARY")
    print("="*80)
    print(f"Visualization approach:")
    print(f"  Species columns: Binary presence (red = present, white = absent)")
    print(f"  Abundance column: Color gradient (white = low, dark blue = high)")
    print(f"  Metaproteome: Split into exact (dark green) and BLAST-only (light green)")
    print(f"  Method: IEDB (purple) vs K-mer (orange)")
    print(f"  Surface accessibility: Accessible (red) vs Not accessible (gray)")
    
    print(f"\nData summary:")
    print(f"  Epitopes: {len(top_epitopes)}")
    print(f"  Abundance range: {top_epitopes['Total_Abundance'].min():.4f} - {top_epitopes['Total_Abundance'].max():.4f}")
    print(f"  Metaproteome exact: {top_epitopes['Metaproteome_Exact'].sum()}")
    print(f"  Metaproteome BLAST: {top_epitopes['Metaproteome_BLAST'].sum()}")
    print(f"  Surface accessible: {top_epitopes['Surface_Accessible'].sum()}")

if __name__ == "__main__":
    metaproteome_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation'
    output_dir = f"{metaproteome_dir}/visualizations"
    log_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs'
    
    create_hybrid_heatmap(metaproteome_dir, output_dir, log_dir)