#!/bin/bash
#SBATCH --job-name=epitope_blast
#SBATCH --account=project_2009813
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/logs/blast_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/logs/blast_%j.err

# Load BLAST module
module load biokit

# Define paths
PROTEIN_DB="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/bEVs_all_sequences.fa"
EPITOPE_QUERY="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/iedb_epitopes.fa"
OUTPUT_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results"
DB_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/blast_db"

# Create output and logs directories if they don't exist
mkdir -p ${OUTPUT_DIR}
mkdir -p /scratch/project_2009813/JOHANNA/II_potential_epitopes/logs
mkdir -p ${DB_DIR}

echo "============================================"
echo "Starting IEDB Epitope BLAST Analysis"
echo "Date: $(date)"
echo "============================================"

# Step 1: Create BLAST database from your proteins
echo "Creating BLAST database..."
makeblastdb -in ${PROTEIN_DB} \
            -dbtype prot \
            -out ${DB_DIR}/bEVs_proteins

echo "Database created successfully!"

# Step 2: Run BLASTP
echo "Running BLASTP..."
blastp -query ${EPITOPE_QUERY} \
       -db ${DB_DIR}/bEVs_proteins \
       -out ${OUTPUT_DIR}/epitope_blast_results.txt \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq" \
       -evalue 0.01 \
       -num_threads 4 \
       -max_target_seqs 10

echo "BLAST search completed!"

# Step 3: Filter for high-quality matches (≥80% identity)
echo "Filtering for high-quality matches (≥80% identity)..."
awk '$3 >= 80' ${OUTPUT_DIR}/epitope_blast_results.txt > ${OUTPUT_DIR}/epitope_blast_high_quality.txt

# Step 4: Count results
TOTAL_HITS=$(wc -l < ${OUTPUT_DIR}/epitope_blast_results.txt)
HIGH_QUAL=$(wc -l < ${OUTPUT_DIR}/epitope_blast_high_quality.txt)

echo "============================================"
echo "BLAST Analysis Complete!"
echo "Total hits: ${TOTAL_HITS}"
echo "High quality hits (≥80% identity): ${HIGH_QUAL}"
echo "============================================"
echo "Results saved to:"
echo "  All hits: ${OUTPUT_DIR}/epitope_blast_results.txt"
echo "  High quality: ${OUTPUT_DIR}/epitope_blast_high_quality.txt"
echo "============================================"
