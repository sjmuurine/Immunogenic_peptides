#!/bin/bash
#SBATCH --job-name=download_human_proteome
#SBATCH --account=project_2009813
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/download_human_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/download_human_%j.err

OUTPUT_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data"

echo "============================================"
echo "Downloading Human Proteome from UniProt"
echo "Date: $(date)"
echo "============================================"

cd ${OUTPUT_DIR}

# Remove any partial downloads
rm -f human_proteome_swissprot.fasta* 

# Download human proteome (Swiss-Prot reviewed entries)
echo "Downloading Swiss-Prot human proteome..."
wget -O human_proteome_swissprot.fasta \
  "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29"

# Check if download was successful
if [ ! -f human_proteome_swissprot.fasta ]; then
    echo "ERROR: Download failed!"
    exit 1
fi

# Check file size
FILE_SIZE=$(stat -f%z human_proteome_swissprot.fasta 2>/dev/null || stat -c%s human_proteome_swissprot.fasta)
echo "File size: ${FILE_SIZE} bytes"

# Count sequences
NUM_SEQS=$(grep -c "^>" human_proteome_swissprot.fasta)

echo "============================================"
echo "Download complete!"
echo "Human proteome sequences: ${NUM_SEQS}"
echo "File: ${OUTPUT_DIR}/human_proteome_swissprot.fasta"
echo "============================================"

# Show first few headers to verify
echo ""
echo "First 5 protein headers:"
grep "^>" human_proteome_swissprot.fasta | head -5