#!/usr/bin/env python3

"""
Corrected Hybrid Heatmap Generator for Epitope Validation
Creates binary species presence heatmap with abundance annotation

Usage: python3 v2_33_create_hybrid_heatmap_fixed.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import os
from datetime import datetime

def create_hybrid_heatmap_fixed():
    """
    Create corrected hybrid heatmap with proper species parsing and improved layout.
    """
    
    print("="*80)
    print("CREATING CORRECTED HYBRID HEATMAP (TOP 60)")
    print("="*80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # File paths for Puhti
    data_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation/REFINED_TIER1_with_metaproteome_validation.txt'
    output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/006_metaproteome_validation/visualizations'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load validated epitope data
    print(f"Loading: {data_file}")
    df = pd.read_csv(data_file, sep='\t')
    print(f"Loaded {len(df)} epitope clusters")
    
    # Get all unique species from the data
    all_species_sets = []
    for _, row in df.iterrows():
        species_str = row.get('Species_List', '')
        if pd.notna(species_str) and species_str != '':
            species = [s.strip() for s in species_str.split(';')]
            all_species_sets.extend(species)
    
    unique_species = sorted(set(all_species_sets))
    print(f"Found species: {unique_species}")
    
    # Define comprehensive species order (Gram+ then Gram-)
    gram_positive = ['B_longum', 'E_faecalis', 'L_casei', 'L_fermentum', 'L_gasseri', 
                    'L_plantarum', 'L_reuteri', 'P_freudenreichii', 'S_salivarius', 
                    'S_parvirubra', 'T_sanguinis']
    
    gram_negative = ['A_muciniphila', 'B_fragilis', 'B_thetaiotaomicron', 'E_cloacae', 
                    'F_nucleatum', 'F_prausnitzii', 'G_adiacens', 'P_russellii', 'V_magna']
    
    # Use only species that appear in our data
    gram_pos_present = [s for s in gram_positive if s in unique_species]
    gram_neg_present = [s for s in gram_negative if s in unique_species]
    species_order = gram_pos_present + gram_neg_present
    
    print(f"Using species order: {species_order}")
    print(f"Gram+: {len(gram_pos_present)}, Gram-: {len(gram_neg_present)}")
    
    # Load surface accessibility data
    try:
        intensity_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/anatomia_protein_intensities.txt'
        intensity_df = pd.read_csv(intensity_file, sep='\t')
        protein_localization_map = dict(zip(intensity_df['Protein_ID'], intensity_df['Localization']))
        print(f"Loaded localization data for {len(protein_localization_map)} proteins")
    except:
        print("Could not load localization data - using dummy values")
        protein_localization_map = {}
    
    # Process epitope data
    heatmap_data = []
    
    for _, epitope in df.iterrows():
        epitope_id = epitope['Cluster_ID']
        sequence = epitope['Representative_Sequence']
        species_str = epitope.get('Species_List', '')
        proteins_str = epitope.get('Source_Proteins_All', '')
        
        # Initialize row data
        row_data = {
            'Epitope_ID': epitope_id,
            'Sequence': sequence,
            'Method': epitope.get('Representative_Method', ''),
            'Human_Identity': epitope.get('Max_Human_Identity', 0),
            'Total_Abundance': float(epitope.get('Total_Relative_Abundance', 0)),
            'Metaproteome_Exact': bool(epitope.get('Found_Exact_Match', False)),
            'Metaproteome_BLAST': bool(epitope.get('Found_BLAST_Match', False)),
            'Surface_Accessible': False
        }
        
        # Initialize species presence (binary: 0 = absent, 1 = present)
        for species in species_order:
            row_data[species] = 0
        
        # Parse species presence from Species_List column (MUCH SIMPLER!)
        if pd.notna(species_str) and species_str != '':
            present_species = [s.strip() for s in species_str.split(';')]
            for species in present_species:
                if species in species_order:
                    row_data[species] = 1
        
        # Determine surface accessibility from proteins
        if pd.notna(proteins_str) and proteins_str != '':
            surface_proteins = 0
            total_anatomia_proteins = 0
            
            for protein_full in proteins_str.split(';'):
                protein_full = protein_full.strip()
                
                # Only check ANATOMIA proteins (WP_) for localization
                if protein_full.startswith('WP_'):
                    protein_id = protein_full.split()[0]
                    localization = protein_localization_map.get(protein_id, '')
                    total_anatomia_proteins += 1
                    
                    if localization in ['Cytoplasmic membrane', 'Outer membrane', 
                                      'Extracellular', 'Cellwall']:
                        surface_proteins += 1
            
            # Determine surface accessibility (majority vote)
            if total_anatomia_proteins > 0:
                row_data['Surface_Accessible'] = (surface_proteins / total_anatomia_proteins) > 0.5
        
        heatmap_data.append(row_data)
    
    # Convert to DataFrame and sort
    heatmap_df = pd.DataFrame(heatmap_data)
    
    # Sort by total abundance (descending) then by species count
    species_counts = heatmap_df[species_order].sum(axis=1)
    heatmap_df['Species_Count'] = species_counts
    heatmap_df = heatmap_df.sort_values(['Total_Abundance', 'Species_Count'], ascending=[False, False])
    
    # Take top 60 (or all if fewer)
    top_epitopes = heatmap_df.head(60).copy()
    
    print(f"\nCreating heatmap for top {len(top_epitopes)} epitopes")
    print(f"Abundance range: {top_epitopes['Total_Abundance'].min():.4f} - {top_epitopes['Total_Abundance'].max():.4f}")
    
    # Create the visualization
    fig = plt.figure(figsize=(20, 24))  # Large for readability
    
    # Calculate layout
    n_species = len(species_order)
    n_annotation_cols = 5  # abundance, meta_exact, meta_blast, method, fish
    
    # Create grid layout
    gs = fig.add_gridspec(1, n_species + n_annotation_cols + 1, 
                         width_ratios=[1]*n_species + [0.8]*n_annotation_cols + [2],  # Last for legend
                         hspace=0, wspace=0.02)
    
    # Main species heatmap
    main_ax = fig.add_subplot(gs[0, :n_species])
    
    # Prepare binary species matrix
    species_matrix = top_epitopes[species_order].values
    epitope_labels = [f"{row['Sequence'][:20]}{'...' if len(row['Sequence']) > 20 else ''}" 
                     for _, row in top_epitopes.iterrows()]
    
    # Plot binary heatmap with clear colors
    # Use white (0) = absent, dark teal (1) = present
    colors = ['white', '#2F5F5F']  # White to dark teal (colorblind-friendly)
    cmap = LinearSegmentedColormap.from_list('binary', colors, N=2)
    
    im = main_ax.imshow(species_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    
    # Add grid lines for clarity
    for i in range(len(epitope_labels)+1):
        main_ax.axhline(y=i-0.5, color='lightgray', linewidth=0.5)
    for i in range(len(species_order)+1):
        main_ax.axvline(x=i-0.5, color='lightgray', linewidth=0.5)
    
    # Species labels
    main_ax.set_xticks(range(len(species_order)))
    main_ax.set_xticklabels([s.replace('_', ' ') for s in species_order], 
                           rotation=45, ha='right', fontsize=11)
    main_ax.set_yticks(range(len(epitope_labels)))
    main_ax.set_yticklabels(epitope_labels, fontsize=9)
    
    # Add Gram+/Gram- separation
    gram_pos_count = len(gram_pos_present)
    if gram_pos_count > 0 and len(gram_neg_present) > 0:
        main_ax.axvline(x=gram_pos_count-0.5, color='black', linewidth=3, alpha=0.8)
    
    main_ax.set_xlim(-0.5, n_species-0.5)
    main_ax.set_ylim(len(epitope_labels)-0.5, -0.5)
    main_ax.set_title('Species Presence', fontsize=12, fontweight='bold', pad=20)
    
    # Add Gram labels
    if gram_pos_count > 0:
        main_ax.text(gram_pos_count/2 - 0.5, -3, 'Gram+', ha='center', va='center', 
                    fontsize=12, fontweight='bold')
    if len(gram_neg_present) > 0:
        main_ax.text(gram_pos_count + len(gram_neg_present)/2 - 0.5, -3, 'Gram-', 
                    ha='center', va='center', fontsize=12, fontweight='bold')
    
    # Annotation columns
    col_names = ['Abundance', 'Meta\nExact', 'Meta\nBLAST', 'Method', 'Fish']
    
    # 1. Abundance column (improved color scaling)
    abundance_ax = fig.add_subplot(gs[0, n_species])
    
    # Get abundance values with safety check
    if 'Total_Abundance' in top_epitopes.columns:
        abundance_values = top_epitopes['Total_Abundance'].values
    else:
        print("Warning: Total_Abundance column not found, using zeros")
        abundance_values = np.zeros(len(top_epitopes))
    
    # Improved abundance scaling with log transformation
    if abundance_values.max() > 0:
        # Log scale for better visualization
        log_abundances = np.log10(abundance_values + 1e-6)
        abundance_normalized = (log_abundances - log_abundances.min()) / (log_abundances.max() - log_abundances.min())
    else:
        abundance_normalized = np.zeros_like(abundance_values)
    
    # Better abundance colormap: white to dark blue
    abundance_cmap = LinearSegmentedColormap.from_list('abundance', 
                                                      ['white', '#deebf7', '#9ecae1', 
                                                       '#4292c6', '#2171b5', '#08519c'], N=100)
    
    for i, norm_val in enumerate(abundance_normalized):
        color = abundance_cmap(norm_val)
        abundance_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, 
                                           edgecolor='lightgray', linewidth=0.5))
    
    abundance_ax.set_xlim(0, 1)
    abundance_ax.set_ylim(0, len(top_epitopes))
    abundance_ax.set_xticks([])
    abundance_ax.set_yticks([])
    abundance_ax.set_title('Abundance', fontsize=10, fontweight='bold', pad=10)
    
    # 2. Metaproteome Exact column
    meta_exact_ax = fig.add_subplot(gs[0, n_species + 1])
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        color = '#228B22' if row['Metaproteome_Exact'] else 'white'  # Forest green
        meta_exact_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, 
                                            edgecolor='lightgray', linewidth=0.5))
    
    meta_exact_ax.set_xlim(0, 1)
    meta_exact_ax.set_ylim(0, len(top_epitopes))
    meta_exact_ax.set_xticks([])
    meta_exact_ax.set_yticks([])
    meta_exact_ax.set_title('Meta\nExact', fontsize=10, fontweight='bold', pad=10)
    
    # 3. Metaproteome BLAST column  
    meta_blast_ax = fig.add_subplot(gs[0, n_species + 2])
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        # Show light green only if BLAST but no exact match
        if row['Metaproteome_BLAST'] and not row['Metaproteome_Exact']:
            color = '#98FB98'  # Light green
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
    method_ax = fig.add_subplot(gs[0, n_species + 3])
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        method = str(row['Method'])
        color = '#800080' if 'IEDB' in method else '#FF8C00'  # Purple vs Dark orange
        method_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color,
                                        edgecolor='lightgray', linewidth=0.5))
    
    method_ax.set_xlim(0, 1)
    method_ax.set_ylim(0, len(top_epitopes))
    method_ax.set_xticks([])
    method_ax.set_yticks([])
    method_ax.set_title('Method', fontsize=10, fontweight='bold', pad=10)
    
    # 5. Surface accessibility (Fish) column
    fish_ax = fig.add_subplot(gs[0, n_species + 4])
    for i, (_, row) in enumerate(top_epitopes.iterrows()):
        color = '#B22222' if row['Surface_Accessible'] else '#696969'  # Firebrick vs Dim gray
        fish_ax.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color,
                                      edgecolor='lightgray', linewidth=0.5))
    
    fish_ax.set_xlim(0, 1)
    fish_ax.set_ylim(0, len(top_epitopes))
    fish_ax.set_xticks([])
    fish_ax.set_yticks([])
    fish_ax.set_title('Fish', fontsize=10, fontweight='bold', pad=10)
    
    # Legend
    legend_ax = fig.add_subplot(gs[0, -1])
    legend_ax.axis('off')
    
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='black', label='Not present'),
        plt.Rectangle((0,0),1,1, facecolor='#2F5F5F', label='Species present'),
        plt.Rectangle((0,0),1,1, facecolor='#08519c', label='High abundance'),
        plt.Rectangle((0,0),1,1, facecolor='white', edgecolor='gray', label='Low abundance'),
        plt.Rectangle((0,0),1,1, facecolor='#228B22', label='Metaproteome: Exact'),
        plt.Rectangle((0,0),1,1, facecolor='#98FB98', label='Metaproteome: BLAST only'),
        plt.Rectangle((0,0),1,1, facecolor='#800080', label='Method: IEDB'),
        plt.Rectangle((0,0),1,1, facecolor='#FF8C00', label='Method: K-mer'),
        plt.Rectangle((0,0),1,1, facecolor='#B22222', label='Surface accessible'),
        plt.Rectangle((0,0),1,1, facecolor='#696969', label='Not surface accessible'),
    ]
    
    legend_ax.legend(handles=legend_elements, loc='center left', fontsize=12,
                    bbox_to_anchor=(0, 0.5))
    
    # Main title
    fig.suptitle(f'TOP{len(top_epitopes)} Epitopes: Multi-Source Validation with Binary Presence + Abundance\n' +
                f'Species Coverage (Gram+ | Gram-)', 
                fontsize=18, fontweight='bold', y=0.95)
    
    # Save figure
    output_file = f"{output_dir}/hybrid_heatmap_corrected_top60.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Saved corrected heatmap: {output_file}")
    plt.close()
    
    # Save data
    data_file = f"{output_dir}/hybrid_heatmap_data_corrected_top60.txt"
    top_epitopes.to_csv(data_file, sep='\t', index=False)
    print(f"✓ Saved data: {data_file}")
    
    # Summary
    print(f"\n" + "="*80)
    print("CORRECTED HEATMAP SUMMARY")
    print("="*80)
    print(f"Epitopes: {len(top_epitopes)}")
    print(f"Species: {len(species_order)} ({len(gram_pos_present)} Gram+, {len(gram_neg_present)} Gram-)")
    print(f"Species coverage per epitope: {species_counts.mean():.1f} ± {species_counts.std():.1f}")
    print(f"Abundance range: {top_epitopes['Total_Abundance'].min():.4f} - {top_epitopes['Total_Abundance'].max():.4f}")
    print(f"Metaproteome exact: {top_epitopes['Metaproteome_Exact'].sum()}")
    print(f"Metaproteome BLAST: {top_epitopes['Metaproteome_BLAST'].sum()}")
    print(f"Surface accessible: {top_epitopes['Surface_Accessible'].sum()}")

if __name__ == "__main__":
    create_hybrid_heatmap_fixed()
