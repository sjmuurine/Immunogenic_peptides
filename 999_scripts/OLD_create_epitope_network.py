import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import os

def create_epitope_network(output_dir):
    """
    Create network diagram with spring layout for natural clustering.
    """
    
    print("="*70)
    print("CREATING EPITOPE NETWORK DIAGRAM")
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
    
    def parse_species_og_from_header(header, data_source):
        if pd.isna(header):
            return None, None
        
        header_str = str(header)
        
        if data_source == 'LITERATURE':
            species_map = {
                'fragilis': 'B_fragilis', 'nucleatum': 'F_nucleatum', 'muciniphila': 'A_muciniphila',
                'prausnitzii': 'F_prausnitzii', 'reuteri': 'L_reuteri', 'cloacae': 'E_cloacae',
                'faecalis': 'E_faecalis', 'salivarius': 'S_salivarius', 'longum': 'B_longum',
                'plantarum': 'L_plantarum', 'casei': 'L_casei', 'sanguinis': 'T_sanguinis', 'magna': 'V_magna'
            }
            species = None
            for key, value in species_map.items():
                if key in header_str.lower():
                    species = value
                    break
            return species, None
        
        if ' ' not in header_str:
            return None, None
        
        parts = header_str.split(' ', 1)[1].split('_')
        og = None
        og_idx = None
        for i, part in enumerate(parts):
            if part.startswith('OG') or part == 'unique':
                og = part
                og_idx = i
                break
        
        species = None
        if og_idx and og_idx > 0:
            species = '_'.join(parts[:og_idx])
        
        return species, og
    
    mapping_df['Parsed_Species'], mapping_df['Parsed_OG'] = zip(*mapping_df.apply(
        lambda row: parse_species_og_from_header(row['Full_Protein_Header'], row['Data_Source']), axis=1
    ))
    
    # Gram classification
    gram_positive = ['B_longum', 'E_faecalis', 'L_casei', 'L_plantarum', 'L_reuteri', 
                     'P_russellii', 'S_salivarius', 'T_sanguinis']
    gram_negative = ['A_muciniphila', 'B_fragilis', 'E_cloacae', 'F_nucleatum', 
                     'F_prausnitzii', 'V_magna']
    
    method_colors = {
        'HIGH_CONFIDENCE': '#E69F00',
        'IEDB': '#56B4E9',
        'Kmer_EXACT': '#009E73',
        'Kmer_FUZZY': '#CC79A7'
    }
    
    def get_method_color(methods_str):
        if pd.isna(methods_str):
            return '#999999'
        m = str(methods_str)
        for key in ['HIGH_CONFIDENCE', 'IEDB', 'Kmer_EXACT']:
            if key in m:
                return method_colors[key]
        return method_colors['Kmer_FUZZY']
    
    # Create network - use TOP 15 for clarity
    print("\nBuilding network...")
    G = nx.Graph()
    
    epitope_nodes = []
    og_nodes = set()
    species_nodes = set()
    
    top_epitopes = epitope_df.sort_values('Priority_Score', ascending=False).head(15)
    
    for _, epitope_row in top_epitopes.iterrows():
        cluster_id = epitope_row['Cluster_ID']
        epitope_seq = epitope_row['Representative_Sequence']
        
        G.add_node(cluster_id, 
                   node_type='epitope',
                   sequence=epitope_seq,
                   method=epitope_row['Methods_Used'],
                   priority=epitope_row['Priority_Score'])
        epitope_nodes.append(cluster_id)
        
        epitope_ids = cluster_to_epitopes.get(cluster_id, [])
        
        for epitope_id in epitope_ids:
            matches = mapping_df[mapping_df['Epitope_ID'] == epitope_id]
            
            for _, match in matches.iterrows():
                species = match['Parsed_Species']
                og = match['Parsed_OG']
                
                if pd.notna(og) and og != 'unique':
                    if og not in G:
                        G.add_node(og, node_type='og')
                    og_nodes.add(og)
                    G.add_edge(cluster_id, og)
                    
                    if pd.notna(species) and species not in ['LITERATURE', 'Unknown']:
                        if species not in G:
                            gram = 'positive' if species in gram_positive else 'negative'
                            G.add_node(species, node_type='species', gram=gram)
                        species_nodes.add(species)
                        G.add_edge(og, species)
    
    print(f"Network: {len(epitope_nodes)} epitopes, {len(og_nodes)} OGs, {len(species_nodes)} species")
    
    # Use spring layout for natural clustering
    print("\nCreating spring layout...")
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(18, 14))
    
    # Draw edges first
    nx.draw_networkx_edges(G, pos, alpha=0.2, width=1, edge_color='gray', ax=ax)
    
    # Draw epitope nodes (triangles) - smaller and uniform size
    epitope_colors = [get_method_color(G.nodes[node].get('method', '')) for node in epitope_nodes]
    
    nx.draw_networkx_nodes(G, pos, nodelist=epitope_nodes,
                          node_shape='^', node_color=epitope_colors,
                          node_size=800, ax=ax, edgecolors='black', linewidths=2)
    
    # Draw OG nodes (circles)
    nx.draw_networkx_nodes(G, pos, nodelist=list(og_nodes),
                          node_shape='o', node_color='#66c2a5',
                          node_size=500, ax=ax, edgecolors='black', linewidths=2)
    
    # Draw species nodes (squares)
    species_list = list(species_nodes)
    species_colors = []
    for node in species_list:
        gram = G.nodes[node].get('gram', 'unknown')
        if gram == 'positive':
            species_colors.append('#8da0cb')
        else:
            species_colors.append('#fc8d62')
    
    nx.draw_networkx_nodes(G, pos, nodelist=species_list,
                          node_shape='s', node_color=species_colors,
                          node_size=600, ax=ax, edgecolors='black', linewidths=2)
    
    # Add labels with offset to avoid overlap
    label_pos = {}
    for node, (x, y) in pos.items():
        node_type = G.nodes[node].get('node_type', 'unknown')
        if node_type == 'epitope':
            label_pos[node] = (x, y - 0.08)  # Below epitope
        elif node_type == 'og':
            label_pos[node] = (x, y + 0.05)  # Above OG
        else:
            label_pos[node] = (x + 0.08, y)  # Right of species
    
    # Create short labels
    labels = {}
    for node in G.nodes():
        node_type = G.nodes[node].get('node_type', 'unknown')
        if node_type == 'epitope':
            seq = G.nodes[node]['sequence']
            labels[node] = seq[:6] + '..' if len(seq) > 6 else seq
        elif node_type == 'og':
            labels[node] = node.replace('OG', '')  # Just the number
        else:
            labels[node] = node.replace('_', '\n')  # Species on two lines
    
    nx.draw_networkx_labels(G, label_pos, labels=labels, font_size=8, 
                           font_weight='bold', ax=ax)
    
    # Title
    ax.set_title('TOP15 Epitopes: Network of Epitope-Orthogroup-Species Relationships',
                fontsize=16, fontweight='bold', pad=20)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=method_colors['HIGH_CONFIDENCE'], label='HIGH_CONFIDENCE', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['IEDB'], label='IEDB', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['Kmer_EXACT'], label='Kmer_EXACT', edgecolor='black'),
        mpatches.Patch(facecolor=method_colors['Kmer_FUZZY'], label='Kmer_FUZZY', edgecolor='black'),
        mpatches.Patch(facecolor='white', label='', edgecolor='none'),
        mpatches.Patch(facecolor='#66c2a5', label='Orthogroup', edgecolor='black'),
        mpatches.Patch(facecolor='white', label='', edgecolor='none'),
        mpatches.Patch(facecolor='#8da0cb', label='Gram+ species', edgecolor='black'),
        mpatches.Patch(facecolor='#fc8d62', label='Gram- species', edgecolor='black'),
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, 
             frameon=True, title='Node Types', title_fontsize=11)
    
    # Add shape legend
    from matplotlib.lines import Line2D
    shape_legend = [
        Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markersize=10, label='Epitope'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Orthogroup'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markersize=10, label='Species'),
    ]
    ax.legend(handles=shape_legend, loc='upper left', fontsize=10,
             frameon=True, title='Node Shapes', title_fontsize=11)
    
    # Add color legend back
    ax2 = ax.twinx()
    ax2.axis('off')
    ax2.legend(handles=legend_elements, loc='lower right', fontsize=10,
              frameon=True, title='Epitope Methods', title_fontsize=11)
    
    ax.axis('off')
    ax.margins(0.1)
    
    # Save
    output_file = f"{output_dir}/network_epitope_og_species_TOP15.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\n✓ Saved: {output_file}")
    plt.close()
    
    print("\n" + "="*70)
    print("NETWORK COMPLETE")
    print("="*70)

# Run
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/007_visualizations'
os.makedirs(output_dir, exist_ok=True)

create_epitope_network(output_dir)