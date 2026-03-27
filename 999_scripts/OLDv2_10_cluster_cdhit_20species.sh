#!/bin/bash
#SBATCH --job-name=cluster_cdhit
#SBATCH --account=project_2009813
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_cdhit_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_cdhit_%j.err

echo "=================================================="
echo "CD-HIT Epitope Clustering (20 species corrected)"
echo "=================================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module load biokit
module load gcc/11.3.0
module load biopythontools

# Directories
BASE_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes"
INPUT_FILE="${BASE_DIR}/003_result_comparisons/epitopes_PASSED_human_filter.txt"
OUTPUT_DIR="${BASE_DIR}/004_cluster"
TEMP_DIR="${OUTPUT_DIR}/temp"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TEMP_DIR}

echo "Step 1: Preparing FASTA file from epitopes"
echo "============================================"

python3 << 'PREP_FASTA'
import pandas as pd

# Load epitopes
input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
output_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/temp/epitopes.fasta'

print(f"Loading: {input_file}")
epitopes_df = pd.read_csv(input_file, sep='\t')

print(f"Total epitopes: {len(epitopes_df):,}")

# Filter by length
epitopes_df = epitopes_df[
    (epitopes_df['Epitope_Sequence'].str.len() >= 11) & 
    (epitopes_df['Epitope_Sequence'].str.len() <= 24)
]

print(f"After length filter (11-24 aa): {len(epitopes_df):,}")

# Create FASTA with index as ID
print(f"Writing FASTA: {output_fasta}")
with open(output_fasta, 'w') as f:
    for idx, row in epitopes_df.iterrows():
        f.write(f">{idx}\n{row['Epitope_Sequence']}\n")

print(f"✓ FASTA created with {len(epitopes_df):,} sequences")

# Save index mapping
mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/temp/index_mapping.txt'
epitopes_df[['Epitope_Sequence', 'Epitope_ID', 'Method', 'Source_Proteins', 'Num_Source_Proteins', 
              'Max_Human_Identity', 'Max_Human_Coverage', 'Data_Source']].to_csv(mapping_file, sep='\t', index=True)
print(f"✓ Index mapping saved: {mapping_file}")

PREP_FASTA

echo ""
echo "Step 2: Running CD-HIT"
echo "======================"

FASTA_FILE="${TEMP_DIR}/epitopes.fasta"
CDHIT_OUTPUT="${TEMP_DIR}/epitopes_cdhit"

echo "Running CD-HIT with 80% identity threshold..."

cd-hit -i ${FASTA_FILE} \
       -o ${CDHIT_OUTPUT} \
       -c 0.8 \
       -n 4 \
       -d 0 \
       -T 4 \
       -M 0 \
       -g 1

if [ $? -eq 0 ]; then
    echo "✓ CD-HIT completed successfully"
else
    echo "✗ CD-HIT failed!"
    exit 1
fi

echo ""
echo "Step 3: Parsing CD-HIT output and aggregating metadata"
echo "======================================================"

python3 << 'PARSE_CLUSTERS'
import pandas as pd
from collections import defaultdict
import re

def parse_species_from_protein_id(protein_id):
    """Parse species from protein ID - ALL 20 SPECIES (CORRECTED)."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    # All 20 species patterns (order matters - check specific before general)
    species_patterns = {
        'thetaiotaomicron': 'B_thetaiotaomicron',
        'fragilis': 'B_fragilis',
        'longum': 'B_longum',
        'nucleatum': 'F_nucleatum',
        'prausnitzii': 'F_prausnitzii',
        'freudenreichii': 'P_freudenreichii',  # Handles P_freudenreichii_s_freudenreichii
        'muciniphila': 'A_muciniphila',
        'reuteri': 'L_reuteri',
        'plantarum': 'L_plantarum',
        'casei': 'L_casei',
        'fermentum': 'L_fermentum',
        'gasseri': 'L_gasseri',
        'cloacae': 'E_cloacae',
        'faecalis': 'E_faecalis',
        'salivarius': 'S_salivarius',
        'sanguinis': 'T_sanguinis',
        'russellii': 'P_russellii',
        'magna': 'V_magna',
        'adiacens': 'G_adiacens',
        'parvirubra': 'S_parvirubra'
    }
    
    # LITERATURE format - case insensitive matching
    if 'LITERATURE' in protein_id:
        protein_lower = protein_id.lower()
        for pattern, species in species_patterns.items():
            if pattern in protein_lower:
                return species
        return None
    
    # ANATOMIA format: "ProteinID Species_OG_Conservation_Localization"
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

# Load index mapping
mapping_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/temp/index_mapping.txt'
print(f"Loading index mapping: {mapping_file}")
mapping_df = pd.read_csv(mapping_file, sep='\t', index_col=0)

# Parse CD-HIT cluster file
cluster_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/temp/epitopes_cdhit.clstr'
print(f"Parsing CD-HIT clusters: {cluster_file}")

clusters = defaultdict(list)
current_cluster = None

with open(cluster_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>Cluster'):
            current_cluster = int(line.split()[1])
        else:
            # Parse sequence ID from cluster line
            match = re.search(r'>(\d+)', line)
            if match:
                seq_idx = int(match.group(1))
                is_representative = '*' in line
                clusters[current_cluster].append((seq_idx, is_representative))

print(f"Found {len(clusters)} clusters")

# Process clusters
cluster_results = []

print("\nProcessing clusters and parsing species (ALL 20 SPECIES)...")
species_found = set()

for cluster_id, members in clusters.items():
    if (cluster_id + 1) % 10000 == 0:
        print(f"  Processing cluster {cluster_id + 1:,}/{len(clusters):,}...", end='\r')
    
    # Get member indices
    member_indices = [m[0] for m in members]
    
    # Find representative (CD-HIT marks it with *)
    rep_idx = [m[0] for m in members if m[1]]
    if rep_idx:
        rep_idx = rep_idx[0]
    else:
        rep_idx = member_indices[0]
    
    # Check if any member has IEDB - if so, prefer that as representative
    iedb_members = [idx for idx in member_indices if 'IEDB' in str(mapping_df.loc[idx, 'Method'])]
    if iedb_members:
        rep_idx = iedb_members[0]
    
    representative = mapping_df.loc[rep_idx]
    
    # Aggregate proteins from all members
    all_proteins = []
    for idx in member_indices:
        proteins_str = str(mapping_df.loc[idx, 'Source_Proteins'])
        if proteins_str and proteins_str != 'nan':
            all_proteins.extend([p.strip() for p in proteins_str.split(';')])
    
    all_proteins = list(set(all_proteins))
    
    # Parse species (CORRECTED PARSER)
    species_set = set()
    for protein in all_proteins:
        species = parse_species_from_protein_id(protein)
        if species:
            species_set.add(species)
            species_found.add(species)
    
    # Aggregate methods
    all_methods = set()
    for idx in member_indices:
        method = str(mapping_df.loc[idx, 'Method'])
        if method and method != 'nan':
            all_methods.add(method)
    
    # Save cluster info
    cluster_results.append({
        'Cluster_ID': f"Cluster_{cluster_id + 1}",
        'Cluster_Size': len(members),
        'Representative_Sequence': representative['Epitope_Sequence'],
        'Representative_ID': representative['Epitope_ID'],
        'Representative_Method': representative['Method'],
        'Methods_All': '; '.join(sorted(all_methods)),
        'Has_IEDB': any('IEDB' in m for m in all_methods),
        'Num_Species': len(species_set),
        'Species_List': '; '.join(sorted(species_set)),
        'Num_Source_Proteins': len(all_proteins),
        'Max_Human_Identity': representative['Max_Human_Identity'],
        'Max_Human_Coverage': representative['Max_Human_Coverage']
    })

print()

# Report species found
print(f"\nSpecies found across all clusters: {len(species_found)}")
print("Species list:")
for species in sorted(species_found):
    print(f"  - {species}")

# Create DataFrame and sort
cluster_df = pd.DataFrame(cluster_results)
cluster_df = cluster_df.sort_values('Num_Species', ascending=False)

# Save
output_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/epitope_clusters_summary.txt'
cluster_df.to_csv(output_file, sep='\t', index=False)

print(f"\n✓ Saved: {output_file}")
print(f"\nCluster statistics:")
print(f"  Total clusters: {len(cluster_df):,}")
print(f"  Singletons: {(cluster_df['Cluster_Size'] == 1).sum():,}")
print(f"  Multi-member: {(cluster_df['Cluster_Size'] > 1).sum():,}")
print(f"  Max species per cluster: {cluster_df['Num_Species'].max()}")
print(f"  Mean species per cluster: {cluster_df['Num_Species'].mean():.1f}")

print(f"\nTop 10 clusters by species breadth:")
for _, row in cluster_df.head(10).iterrows():
    print(f"  {row['Representative_Sequence'][:40]:40} - {row['Num_Species']:2} species ({row['Species_List'][:60]}...)")

PARSE_CLUSTERS

echo ""
echo "Step 4: Cleaning up temporary files"
echo "===================================="

# Keep cluster files for reference, remove only FASTA
rm ${TEMP_DIR}/epitopes.fasta

echo "✓ Cleanup complete"
echo ""
echo "=================================================="
echo "Clustering Complete!"
echo "End time: $(date)"
echo "=================================================="