#!/bin/bash
#SBATCH --job-name=epitope_blast_both
#SBATCH --account=project_2009813
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/blast_both_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/blast_both_%j.err

# Load BLAST module
module load biokit

# Define paths
ANATOMIA_DB="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/ANATOMIA_bEVs_all_sequences.fa"
LITERATURE_DB="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/literature_merged.fasta"
EPITOPE_QUERY="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/iedb_epitopes.fa"
OUTPUT_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results"
DB_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/blast_db"

# Create directories if they don't exist
mkdir -p ${DB_DIR}

echo "============================================"
echo "IEDB Epitope BLAST: ANATOMIA + Literature"
echo "Date: $(date)"
echo "============================================"

# Step 1: Create BLAST databases
echo ""
echo "Creating BLAST database for ANATOMIA data..."
makeblastdb -in ${ANATOMIA_DB} \
            -dbtype prot \
            -out ${DB_DIR}/anatomia_db

echo "Creating BLAST database for Literature data..."
makeblastdb -in ${LITERATURE_DB} \
            -dbtype prot \
            -out ${DB_DIR}/literature_db

echo "Databases created successfully!"

# Step 2: BLAST against ANATOMIA data
echo ""
echo "BLASTing epitopes against ANATOMIA data..."
blastp -query ${EPITOPE_QUERY} \
       -db ${DB_DIR}/anatomia_db \
       -out ${OUTPUT_DIR}/epitope_blast_ANATOMIA.txt \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq" \
       -evalue 0.01 \
       -num_threads 4 \
       -max_target_seqs 10

echo "ANATOMIA BLAST completed!"

# Step 3: BLAST against Literature data
echo ""
echo "BLASTing epitopes against Literature data..."
blastp -query ${EPITOPE_QUERY} \
       -db ${DB_DIR}/literature_db \
       -out ${OUTPUT_DIR}/epitope_blast_LITERATURE.txt \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq" \
       -evalue 0.01 \
       -num_threads 4 \
       -max_target_seqs 10

echo "Literature BLAST completed!"

# Step 4: Add source labels and combine results
echo ""
echo "Combining results..."

# Add "ANATOMIA" label to first dataset
awk '{print $0"\tANATOMIA"}' ${OUTPUT_DIR}/epitope_blast_ANATOMIA.txt > ${OUTPUT_DIR}/temp_anatomia.txt

# Add "LITERATURE" label to second dataset
awk '{print $0"\tLITERATURE"}' ${OUTPUT_DIR}/epitope_blast_LITERATURE.txt > ${OUTPUT_DIR}/temp_literature.txt

# Combine both files with header
echo -e "Epitope_ID\tProtein_ID\tPercent_Identity\tAlignment_Length\tEpitope_Length\tProtein_Length\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tEpitope_Sequence\tProtein_Sequence\tData_Source" > ${OUTPUT_DIR}/epitope_blast_COMBINED.txt
cat ${OUTPUT_DIR}/temp_anatomia.txt ${OUTPUT_DIR}/temp_literature.txt >> ${OUTPUT_DIR}/epitope_blast_COMBINED.txt

# Clean up temp files
rm ${OUTPUT_DIR}/temp_anatomia.txt ${OUTPUT_DIR}/temp_literature.txt

echo "Combined results created!"

# Step 5: Create filtered high-quality results (≥80% identity)
echo ""
echo "Filtering for high-quality matches (≥80% identity)..."

# Filter combined results
awk 'NR==1 || $3 >= 80' ${OUTPUT_DIR}/epitope_blast_COMBINED.txt > ${OUTPUT_DIR}/epitope_blast_COMBINED_high_quality.txt

# Also create separate high-quality files
awk '$3 >= 80' ${OUTPUT_DIR}/epitope_blast_ANATOMIA.txt > ${OUTPUT_DIR}/epitope_blast_ANATOMIA_high_quality.txt
awk '$3 >= 80' ${OUTPUT_DIR}/epitope_blast_LITERATURE.txt > ${OUTPUT_DIR}/epitope_blast_LITERATURE_high_quality.txt

# Step 6: Generate summary statistics
echo ""
echo "============================================"
echo "BLAST RESULTS SUMMARY"
echo "============================================"

ANATOMIA_TOTAL=$(wc -l < ${OUTPUT_DIR}/epitope_blast_ANATOMIA.txt)
ANATOMIA_HQ=$(wc -l < ${OUTPUT_DIR}/epitope_blast_ANATOMIA_high_quality.txt)

LITERATURE_TOTAL=$(wc -l < ${OUTPUT_DIR}/epitope_blast_LITERATURE.txt)
LITERATURE_HQ=$(wc -l < ${OUTPUT_DIR}/epitope_blast_LITERATURE_high_quality.txt)

COMBINED_TOTAL=$((ANATOMIA_TOTAL + LITERATURE_TOTAL))
COMBINED_HQ=$((ANATOMIA_HQ + LITERATURE_HQ))

echo ""
echo "ANATOMIA Dataset:"
echo "  Total hits: ${ANATOMIA_TOTAL}"
echo "  High quality (≥80%): ${ANATOMIA_HQ}"
echo ""
echo "Literature Dataset:"
echo "  Total hits: ${LITERATURE_TOTAL}"
echo "  High quality (≥80%): ${LITERATURE_HQ}"
echo ""
echo "Combined:"
echo "  Total hits: ${COMBINED_TOTAL}"
echo "  High quality (≥80%): ${COMBINED_HQ}"
echo ""
echo "============================================"
echo "Output files:"
echo "  Combined results: epitope_blast_COMBINED.txt"
echo "  Combined HQ: epitope_blast_COMBINED_high_quality.txt"
echo "  ANATOMIA results: epitope_blast_ANATOMIA.txt"
echo "  ANATOMIA HQ: epitope_blast_ANATOMIA_high_quality.txt"
echo "  Literature results: epitope_blast_LITERATURE.txt"
echo "  Literature HQ: epitope_blast_LITERATURE_high_quality.txt"
echo "============================================"
echo "Analysis complete: $(date)"
echo "============================================"
