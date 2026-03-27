#!/bin/bash
#SBATCH --job-name=cluster_cdhit_safe
#SBATCH --account=project_2009813
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_cdhit_safe_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/cluster_cdhit_safe_%j.err

echo "=================================================="
echo "CD-HIT Clustering with SAFE Representative Selection"
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

input_file = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons/epitopes_PASSED_human_filter.txt'
output_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster/temp/epitopes.fasta'

print(f"Loading: {input_file}")
epitopes_df = pd.read_csv(input_file, sep='\t')
print(f"Total epitopes: {len(epitopes_df):,}")

epitopes_df = epitopes_df[
    (epitopes_df['Epitope_Sequence'].str.len() >= 11) & 
    (epitopes_df['Epitope_Sequence'].str.len() <= 24)
]
print(f"After length filter (11-24 aa): {len(epitopes_df):,}")

print(f"Writing FASTA: {output_fasta}")
with open(output_fasta, 'w') as f:
    for idx, row in epitopes_df.iterrows():
        f.write(f">{idx}\n{row['Epitope_Sequence']}\n")

print(f"✓ FASTA created with {len(epitopes_df):,} sequences")

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
echo "Step 3: SAFE Representative Selection + Species Aggregation"
echo "=========================================================="

python3 << 'PARSE_CLUSTERS'
import pandas as pd
from collections import defaultdict
import re

def parse_species_from_protein_id(protein_id):
    """Parse species from protein ID - 20 SPECIES (B_longum_suis merged to B_longum)."""
    if pd.isna(protein_id):
        return None
    
    protein_id = str(protein_id)
    
    species_patterns = {
        'thetaiotaomicron': 'B_thetaiotaomicron',
        'fragilis': 'B_fragilis',
        'longum': 'B_longum',  # Handles both B_longum and B_longum_suis
        'nucleatum': 'F_nucleatum',
        'prausnitzii': 'F_prausnitzii',
        'freudenreichii': 'P_freudenreichii',
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
    
    if 'LITERATURE' in protein_id:
        protein_lower = protein_id.lower()
        for pattern, species in species_patterns.items():
            if pattern in protein_lower:
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
        species = '_'.join(parts[:og_idx])
        if species == 'B_longum_suis':
            species = 'B_longum'
        return species
    
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
            match = re.search(r'>(\d+)', line)
            if match:
                seq_idx = int(match.group(1))
                is_representative = '*' in line
                clusters[current_cluster].append((seq_idx, is_representative))

print(f"Found {len(clusters)} clusters")

# Process clusters with SAFE representative selection
cluster_results = []
print("\nProcessing clusters with SAFE representative selection...")

for cluster_id, members in clusters.items():
    if (cluster_id + 1) % 10000 == 0:
        print(f"  Processing cluster {cluster_id + 1:,}/{len(clusters):,}...", end='\r')
    
    member_indices = [m[0] for m in members]
    
    # SAFE REPRESENTATIVE SELECTION - Priority order:
    # 1. Lowest human identity (SAFETY FIRST!)
    # 2. Has IEDB (bonus, but not required)
    # 3. Most source proteins (tiebreaker)
    
    best_idx = None
    best_identity = float('inf')
    best_has_iedb = False
    best_proteins = 0
    
    for idx in member_indices:
        epitope = mapping_df.loc[idx]
        
        identity = epitope.get('Max_Human_Identity', 100)  # Default to 100 if missing
        has_iedb = 'IEDB' in str(epitope.get('Method', ''))
        num_proteins = epitope.get('Num_Source_Proteins', 0)
        
        # Selection logic: prioritize lowest identity, then IEDB, then proteins
        is_better = False
        
        if identity < best_identity:
            is_better = True
        elif identity == best_identity:
            if has_iedb and not best_has_iedb:
                is_better = True
            elif has_iedb == best_has_iedb and num_proteins > best_proteins:
                is_better = True
        
        if is_better:
            best_idx = idx
            best_identity = identity
            best_has_iedb = has_iedb
            best_proteins = num_proteins
    
    representative = mapping_df.loc[best_idx]
    
    # Aggregate proteins from all cluster members
    all_proteins = []
    for idx in member_indices:
        proteins_str = str(mapping_df.loc[idx, 'Source_Proteins'])
        if proteins_str and proteins_str != 'nan':
            all_proteins.extend([p.strip() for p in proteins_str.split(';')])
    
    all_proteins = list(set(all_proteins))
    
    # Parse species
    species_set = set()
    for protein in all_proteins:
        species = parse_species_from_protein_id(protein)
        if species:
            species_set.add(species)
    
    # Aggregate methods
    all_methods = set()
    for idx in member_indices:
        method = str(mapping_df.loc[idx, 'Method'])
        if method and method != 'nan':
            all_methods.add(method)
    
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
print(f"  Max species: {cluster_df['Num_Species'].max()}")
print(f"  Mean species: {cluster_df['Num_Species'].mean():.1f}")

print(f"\nHuman similarity in representatives (SHOULD BE LOWER NOW!):")
print(f"  Mean identity: {cluster_df['Max_Human_Identity'].mean():.1f}%")
print(f"  Median identity: {cluster_df['Max_Human_Identity'].median():.1f}%")
print(f"  Max identity: {cluster_df['Max_Human_Identity'].max():.1f}%")
print(f"  Min identity: {cluster_df['Max_Human_Identity'].min():.1f}%")

# Show improvement
high_identity = (cluster_df['Max_Human_Identity'] >= 95).sum()
print(f"  Representatives with ≥95% identity: {high_identity} (should be much lower!)")

print(f"\nTop 10 clusters by species breadth (with safer representatives):")
for _, row in cluster_df.head(10).iterrows():
    print(f"  {row['Representative_Sequence'][:30]:30} - {row['Num_Species']:2} sp, {row['Max_Human_Identity']:5.1f}% ID")

PARSE_CLUSTERS

echo ""
echo "Step 4: Cleaning up"
echo "==================="
rm ${TEMP_DIR}/epitopes.fasta
echo "✓ Cleanup complete"

echo ""
echo "=================================================="
echo "CLUSTERING WITH SAFE REPRESENTATIVES COMPLETE!"
echo "End time: $(date)"
echo "=================================================="