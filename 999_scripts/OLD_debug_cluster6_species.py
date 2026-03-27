import pandas as pd

print("="*70)
print("INVESTIGATING TPFFNGYRPQFYFRT SPECIES")
print("="*70)

# 1. Check what TOP50_FINAL says
top50 = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_FINAL_with_metadata.txt', sep='\t')
cluster6 = top50[top50['Cluster_ID'] == 'Cluster_6']

print("\nCluster_6 (TPFFNGYRPQFYFRT) in TOP50_FINAL_with_metadata.txt:")
print(f"  Num_Species: {cluster6['Num_Species'].values[0]}")
print(f"  Species_List: {cluster6['Species_List'].values[0]}")

# 2. Check the protein mapping file being used for heatmap
mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/TOP50_protein_mapping.txt'
mapping_df = pd.read_csv(mapping_file, sep='\t')

print(f"\n\nProtein mapping file date: Jan 13 (OLD parsing)")
print(f"Total records: {len(mapping_df)}")
print(f"Columns: {mapping_df.columns.tolist()}")

# 3. Check what species are in OLD mapping for Cluster_6
detailed = pd.read_csv('/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/TOP50_cluster_detailed_v2.txt', sep='\t')
cluster6_epitopes = detailed[detailed['Cluster_ID'] == 'Cluster_6']['Epitope_ID'].tolist()

print(f"\nCluster_6 contains epitope IDs: {cluster6_epitopes}")

cluster6_proteins = mapping_df[mapping_df['Epitope_ID'].isin(cluster6_epitopes)]
print(f"\nProteins for Cluster_6 in OLD mapping: {len(cluster6_proteins)}")

# Check parsed species
if 'Parsed_Species' in cluster6_proteins.columns:
    old_species = cluster6_proteins['Parsed_Species'].dropna().unique()
    print(f"\nSpecies in OLD mapping (Parsed_Species column): {len(old_species)}")
    print(f"  {sorted(old_species)}")

# 4. NOW reparse with CORRECTED function
def parse_species_from_protein_id(protein_id):
    """CORRECTED VERSION."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    species_patterns = {
        'fragilis': 'B_fragilis', 'nucleatum': 'F_nucleatum', 'muciniphila': 'A_muciniphila',
        'prausnitzii': 'F_prausnitzii', 'reuteri': 'L_reuteri', 'cloacae': 'E_cloacae',
        'faecalis': 'E_faecalis', 'salivarius': 'S_salivarius', 'longum': 'B_longum',
        'plantarum': 'L_plantarum', 'casei': 'L_casei', 'sanguinis': 'T_sanguinis',
        'magna': 'V_magna', 'russellii': 'P_russellii', 'fermentum': 'L_fermentum',
        'gasseri': 'L_gasseri', 'freudenreichii': 'P_freudenreichii'
    }
    
    if 'LITERATURE' in protein_id:
        for pattern, species in species_patterns.items():
            if pattern.lower() in protein_id.lower():
                return species
        return None
    
    if ' ' not in protein_id:
        return None
    
    parts = protein_id.split(' ', 1)[1].split('_')
    og_idx = None
    for i, part in enumerate(parts):
        if part.startswith('OG') or part == 'unique':
            og_idx = i
            break
    
    if og_idx and og_idx > 0:
        return '_'.join(parts[:og_idx])
    
    return None

print("\n" + "="*70)
print("REPARSING WITH CORRECTED FUNCTION")
print("="*70)

cluster6_proteins['Species_NEW'] = cluster6_proteins['Full_Protein_Header'].apply(parse_species_from_protein_id)

new_species = cluster6_proteins['Species_NEW'].dropna().unique()
print(f"\nSpecies with NEW parsing: {len(new_species)}")
print(f"  {sorted(new_species)}")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

print("\nThe heatmap uses the OLD protein mapping file from Jan 13.")
print("This file has the OLD species parsing (before we fixed LITERATURE support).")
print("\nWe need to:")
print("  1. Regenerate TOP50_protein_mapping.txt with corrected parsing")
print("  2. Then regenerate the heatmap")
print("\nOR just update the heatmap script to reparse species on-the-fly!")

print("="*70)