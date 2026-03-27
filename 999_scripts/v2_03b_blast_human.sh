#!/bin/bash
#SBATCH --job-name=human_blast_v2
#SBATCH --account=project_2009813
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/human_blast_v2_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/human_blast_v2_%j.err

echo "=================================================="
echo "Human Proteome BLAST - All Epitopes"
echo "=================================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module load biokit

# Directories
BASE_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes"
WORK_DIR="${BASE_DIR}/003_result_comparisons"
BLAST_DB="${BASE_DIR}/001_data/blast_db/human_proteome_db"

# Input files
INPUT_FASTA="${WORK_DIR}/all_epitopes_for_blast.fasta"
ID_MAPPING="${WORK_DIR}/epitope_id_mapping.txt"

# Output files
BLAST_OUTPUT="${WORK_DIR}/epitopes_vs_human_raw.txt"
BLAST_HEADER="${WORK_DIR}/epitopes_vs_human_with_header.txt"
FINAL_OUTPUT="${WORK_DIR}/epitopes_vs_human_ALL_hits.txt"

echo "Input files:"
echo "  FASTA: ${INPUT_FASTA}"
echo "  ID mapping: ${ID_MAPPING}"
echo "  BLAST DB: ${BLAST_DB}"
echo ""

# Check if input exists
if [ ! -f "${INPUT_FASTA}" ]; then
    echo "ERROR: Input FASTA not found: ${INPUT_FASTA}"
    echo "Please run v2_03a_collect_all_epitopes.py first!"
    exit 1
fi

# Count sequences
NUM_EPITOPES=$(grep -c "^>" "${INPUT_FASTA}")
echo "Total epitopes to BLAST: ${NUM_EPITOPES}"
echo ""

# Run BLAST
echo "Running BLAST against human proteome..."
echo "  Parameters:"
echo "    - E-value threshold: 1000 (permissive)"
echo "    - Max target seqs: 100"
echo "    - Threads: 8"
echo ""

blastp \
    -query "${INPUT_FASTA}" \
    -db "${BLAST_DB}" \
    -out "${BLAST_OUTPUT}" \
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
    -num_threads 8 \
    -evalue 1000 \
    -max_target_seqs 100

# Check if BLAST succeeded
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: BLAST failed!"
    exit 1
fi

echo ""
echo "BLAST completed successfully!"

# Count hits
NUM_HITS=$(wc -l < "${BLAST_OUTPUT}")
echo "  Total hits: ${NUM_HITS}"
echo ""

# Add header to BLAST output
echo "Adding header to BLAST results..."
echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tQuery_Length\tSubject_Length\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tQuery_Coverage\tQuery_Sequence\tSubject_Sequence" > "${BLAST_HEADER}"
cat "${BLAST_OUTPUT}" >> "${BLAST_HEADER}"

echo "  Created: ${BLAST_HEADER}"
echo ""

# Map Query IDs back to epitope sequences
echo "Mapping Query IDs to epitope sequences..."

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import sys

# Paths
mapping_file = sys.argv[1]
blast_file = sys.argv[2]
output_file = sys.argv[3]

print(f"  Loading ID mapping: {mapping_file}")
id_map = pd.read_csv(mapping_file, sep='\t')
id_to_seq = dict(zip(id_map['ID'].astype(str), id_map['Epitope_Sequence']))

print(f"  Loading BLAST results: {blast_file}")
blast_df = pd.read_csv(blast_file, sep='\t')

print(f"  Mapping {len(blast_df)} BLAST hits...")

# Map Query_ID to epitope sequence
blast_df['Epitope_Sequence'] = blast_df['Query_ID'].astype(str).map(id_to_seq)

# Check for unmapped
unmapped = blast_df['Epitope_Sequence'].isna().sum()
if unmapped > 0:
    print(f"  WARNING: {unmapped} hits could not be mapped!")
else:
    print(f"  All hits successfully mapped!")

# Reorder columns
cols = ['Epitope_Sequence', 'Subject_ID', 'Percent_Identity', 'Alignment_Length', 
        'Query_Length', 'Subject_Length', 'Query_Start', 'Query_End', 
        'Subject_Start', 'Subject_End', 'E_value', 'Bit_Score', 'Query_Coverage',
        'Query_Sequence', 'Subject_Sequence']
blast_df = blast_df[cols]

# Save
blast_df.to_csv(output_file, sep='\t', index=False)
print(f"  Saved: {output_file}")
print(f"  Total rows: {len(blast_df)}")

# Summary statistics
print(f"\n  Summary:")
print(f"    Unique epitopes with hits: {blast_df['Epitope_Sequence'].nunique()}")
print(f"    Average hits per epitope: {len(blast_df) / blast_df['Epitope_Sequence'].nunique():.1f}")
print(f"    Max identity: {blast_df['Percent_Identity'].max():.1f}%")
print(f"    Max coverage: {blast_df['Query_Coverage'].max():.1f}%")

PYTHON_SCRIPT

python3 -c "
import sys
sys.argv = ['', '${ID_MAPPING}', '${BLAST_HEADER}', '${FINAL_OUTPUT}']
exec(open('/dev/stdin').read())
" << 'PYTHON_SCRIPT'
import pandas as pd
import sys

mapping_file = sys.argv[1]
blast_file = sys.argv[2]
output_file = sys.argv[3]

print(f"  Loading ID mapping: {mapping_file}")
id_map = pd.read_csv(mapping_file, sep='\t')
id_to_seq = dict(zip(id_map['ID'].astype(str), id_map['Epitope_Sequence']))

print(f"  Loading BLAST results: {blast_file}")
blast_df = pd.read_csv(blast_file, sep='\t')

print(f"  Mapping {len(blast_df)} BLAST hits...")
blast_df['Epitope_Sequence'] = blast_df['Query_ID'].astype(str).map(id_to_seq)

unmapped = blast_df['Epitope_Sequence'].isna().sum()
if unmapped > 0:
    print(f"  WARNING: {unmapped} hits could not be mapped!")
else:
    print(f"  All hits successfully mapped!")

cols = ['Epitope_Sequence', 'Subject_ID', 'Percent_Identity', 'Alignment_Length', 
        'Query_Length', 'Subject_Length', 'Query_Start', 'Query_End', 
        'Subject_Start', 'Subject_End', 'E_value', 'Bit_Score', 'Query_Coverage',
        'Query_Sequence', 'Subject_Sequence']
blast_df = blast_df[cols]

blast_df.to_csv(output_file, sep='\t', index=False)
print(f"  Saved: {output_file}")
print(f"  Total rows: {len(blast_df)}")

print(f"\n  Summary:")
print(f"    Unique epitopes with hits: {blast_df['Epitope_Sequence'].nunique()}")
print(f"    Average hits per epitope: {len(blast_df) / blast_df['Epitope_Sequence'].nunique():.1f}")
print(f"    Max identity: {blast_df['Percent_Identity'].max():.1f}%")
print(f"    Max coverage: {blast_df['Query_Coverage'].max():.1f}%")
PYTHON_SCRIPT

echo ""
echo "Cleaning up temporary files..."
rm "${BLAST_OUTPUT}"
echo "  Removed: ${BLAST_OUTPUT}"

echo ""
echo "=================================================="
echo "BLAST Complete!"
echo "End time: $(date)"
echo "=================================================="
echo ""
echo "Output file: ${FINAL_OUTPUT}"
echo ""
echo "Next step: Run v2_03b_filter_human_similarity.py"
echo "=================================================="