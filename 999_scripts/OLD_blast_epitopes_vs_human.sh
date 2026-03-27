#!/bin/bash
#SBATCH --job-name=epitopes_vs_human
#SBATCH --account=project_2009813
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/human_blast_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/human_blast_%j.err

# Load modules (gcc first!)
module load gcc/11.3.0
module load biokit
module load biopythontools

# Define paths
HUMAN_PROTEOME="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/human_proteome_swissprot.fasta"
COMPARISON_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons"
OUTPUT_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons"
TEMP_DIR="${OUTPUT_DIR}/temp"
DB_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/blast_db"
SCRIPT_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/999_scripts"

mkdir -p ${TEMP_DIR}
mkdir -p ${DB_DIR}

echo "============================================"
echo "BLAST Epitopes vs Human Proteome"
echo "Date: $(date)"
echo "============================================"

# Step 1: Create combined epitope FASTA
echo ""
echo "Step 1: Creating combined epitope FASTA..."
python3 ${SCRIPT_DIR}/create_combined_epitope_fasta.py

if [ ! -f "${TEMP_DIR}/all_epitopes_combined.fasta" ]; then
    echo "ERROR: Failed to create epitope FASTA!"
    exit 1
fi

echo "✓ Combined FASTA created"

# Step 2: Create BLAST database from human proteome
echo ""
echo "Step 2: Creating BLAST database from human proteome..."
makeblastdb -in ${HUMAN_PROTEOME} \
            -dbtype prot \
            -out ${DB_DIR}/human_proteome_db

echo "✓ Database created"

# Step 3: BLAST epitopes against human proteome
echo ""
echo "Step 3: BLASTing epitopes against human proteome..."
blastp -query ${TEMP_DIR}/all_epitopes_combined.fasta \
       -db ${DB_DIR}/human_proteome_db \
       -out ${OUTPUT_DIR}/epitopes_vs_human_ALL_hits.txt \
       -outfmt "6 qseqid stitle pident length qlen slen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
       -num_threads 8 \
       -max_target_seqs 5

echo "✓ BLAST completed"

# Step 4: Filter for bacteria-specific epitopes
echo ""
echo "Step 4: Filtering for bacteria-specific epitopes..."
python3 ${SCRIPT_DIR}/filter_human_blast.py

echo ""
echo "============================================"
echo "ANALYSIS COMPLETE"
echo "Date: $(date)"
echo "============================================"
echo "Results saved to: ${OUTPUT_DIR}"
echo "============================================"