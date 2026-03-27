# Immunogenic_peptides

This is a computational workflow for identifying and validating bacterial epitopes from gut microbiome proteome data with rigorous human cross-reactivity screening.

## Overview

This pipeline processes bacterial protein data to identify potential immunogenic epitopes suitable for research applications. The workflow includes multi-source epitope extraction, human proteome similarity filtering, clustering, abundance analysis, and metaproteome validation.

## Key Features

- **Multi-source epitope identification**: Combines k-mer analysis and IEDB database mining
- **Human cross-reactivity screening**: BLAST-based filtering to avoid autoimmune potential
- **Species-aware clustering**: CD-HIT clustering with representative sequence selection
- **Protein abundance integration**: ANATOMIA (=our lab) proteome intensity data incorporation
- **Metaproteome validation**: Real-world detection confirmation using fecal metaproteome data
- **Multi-criteria prioritization**: Balanced scoring considering species breadth, human similarity, and abundance

## Pipeline Workflow

### Phase 1: Data Collection & Epitope Extraction
1. **K-mer epitope generation** (`999_scripts/v2_01_*.py`)
   - Exact and fuzzy k-mer extraction from bacterial proteins
   - Length filtering (11-24 amino acids)
   
2. **IEDB integration** (`999_scripts/v2_02_*.py`)
   - Incorporation of experimentally validated epitopes
   - Cross-reference with generated k-mers

### Phase 2: Human Cross-Reactivity Screening
3. **Human proteome BLAST** (`999_scripts/v2_03b_blast_human.sh`)
   - Comprehensive similarity screening against human proteome
   - Identity and coverage threshold filtering (<80% identity)
   
4. **Filtering** (`999_scripts/v2_17_human_filter_80percent.py`)
   - Binary filtering to remove potentially cross-reactive sequences
   - Balanced approach between safety and candidate retention

### Phase 3: Clustering & Representative Selection
5. **CD-HIT clustering** (`999_scripts/v2_18_cluster_80pct.sh`)
   - 80% identity clustering to reduce redundancy
   - Selection (lowest human similarity)
   
6. **Low-complexity filtering** (`999_scripts/v2_19_filter_low_complexity.py`)
   - Removal of repetitive/low-complexity sequences
   - Quality control for meaningful epitopes

### Phase 4: Multi-Source Validation & Prioritization
7. **Abundance integration** (`999_scripts/v2_24_merge_intensity_fixed.py`)
   - Integration of ANATOMIA protein abundance data
   - Total relative abundance calculations
   
8. **Metaproteome validation** (`999_scripts/v2_30_metaproteome_validation.py`)
   - BLAST and exact string matching against fecal metaproteome
   - Real-world detection confirmation
   
9. **Multi-criteria prioritization** (`999_scripts/v2_22_prioritize_multi_species.py`)
   - Species breadth and human similarity scoring
   - TIER1/2 candidate selection

### Phase 5: Visualization & Analysis
10. **Comprehensive visualization** (`999_scripts/v2_33_create_hybrid_heatmap_fixed.py`)
    - Multi-source validation heatmaps
    - Species presence, abundance, and validation status

## Directory Structure
```
├── 000_logs/                          # Processing logs
├── 001_data/                          # Input datasets (ANATOMIA, metaproteome)
├── 002_results/                       # Initial analysis results
├── 003_result_comparisons/             # Human similarity analysis results
├── 004_cluster/                       # Clustering results and representatives
├── 005_intensity_analysis/            # Protein abundance integration
├── 006_metaproteome_validation/       # Metaproteome validation results
├── 007_visualizations/                # Final plots and heatmaps
└── 999_scripts/                       # All pipeline scripts
```

## Key Results

- Processed 218,000+ initial epitopes from 15-20 gut bacterial species
- Applied stringent human cross-reactivity filtering (<80% identity)
- Generated refined TIER1 candidates with multi-source validation
- Integrated protein abundance and real-world metaproteome detection
- Produced publication-ready visualizations and candidate prioritization

## Dependencies

- Python 3.6+
- BioPython
- pandas, numpy, matplotlib, seaborn
- CD-HIT
- BLAST+ suite
- SLURM (for HPC execution)

## Usage

Scripts are designed for SLURM-based HPC environments. Key entry points:
```bash
# Human similarity filtering
python3 999_scripts/v2_17_human_filter_80percent.py

# Clustering with safe representatives
sbatch 999_scripts/v2_18_cluster_80pct.sh

# Final prioritization and visualization
python3 999_scripts/v2_22_prioritize_multi_species.py
python3 999_scripts/v2_33_create_hybrid_heatmap_fixed.py
```

## Citation

If you use this pipeline in your research, please cite:
[Publication details to be added upon acceptance]

## Contact

johanna.muurinen@onehealth.fi

---

**Note**: Large data files (>50MB) are not included in this repository due to GitHub size limitations. The pipeline structure and all analysis scripts are provided for reproducibility.
