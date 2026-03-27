#!/bin/bash
# add_headers.sh

RESULTS_DIR="/scratch/project_2009813/JOHANNA/II_potential_epitopes/002_results"

HEADER="Epitope_ID\tProtein_ID\tPercent_Identity\tAlignment_Length\tEpitope_Length\tProtein_Length\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tEpitope_Sequence\tProtein_Sequence"

# Add headers
echo -e "${HEADER}" | cat - ${RESULTS_DIR}/epitope_blast_results.txt > ${RESULTS_DIR}/epitope_blast_results_with_header.txt
echo -e "${HEADER}" | cat - ${RESULTS_DIR}/epitope_blast_high_quality.txt > ${RESULTS_DIR}/epitope_blast_high_quality_with_header.txt

echo "Headers added successfully!"