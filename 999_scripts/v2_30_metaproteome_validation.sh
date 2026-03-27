#!/bin/bash
#SBATCH --job-name=metaproteome_validation
#SBATCH --account=project_2009813
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --output=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/metaproteome_validation_%j.out
#SBATCH --error=/scratch/project_2009813/JOHANNA/II_potential_epitopes/000_logs/metaproteome_validation_%j.err

echo "=================================================="
echo "Metaproteome Validation Analysis"
echo "=================================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""

# Load modules
module load biokit
module load gcc/11.3.0
module load biopythontools

# Run validation
python3 /scratch/project_2009813/JOHANNA/II_potential_epitopes/999_scripts/v2_30_metaproteome_validation.py

echo ""
echo "=================================================="
echo "Metaproteome Validation Complete"
echo "End time: $(date)"
echo "=================================================="