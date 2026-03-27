#!/usr/bin/env python3

def split_fasta_by_og(input_fasta, output_directory):
    """
    Splits a FASTA file into separate files per orthogroup.
    Assumes header format: >OG_ConservedIn_FastaHeaders_Localization
    Skips sequences where OG = "unique"
    """
    
    import os
    
    # Create output directory
    os.makedirs(output_directory, exist_ok=True)
    
    # Dictionary to hold sequences by OG
    og_sequences = {}
    
    # Read the big FASTA
    with open(input_fasta, 'r') as f:
        current_header = None
        current_sequence = []
        
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header and current_sequence:
                    # Extract OG from header (first part before first underscore)
                    og = current_header.split('_')[0].replace('>', '')
                    
                    # Skip if "unique"
                    if og.lower() != "unique":
                        if og not in og_sequences:
                            og_sequences[og] = []
                        og_sequences[og].append((current_header, ''.join(current_sequence)))
                
                # Start new sequence
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Don't forget the last sequence
        if current_header and current_sequence:
            og = current_header.split('_')[0].replace('>', '')
            if og.lower() != "unique":
                if og not in og_sequences:
                    og_sequences[og] = []
                og_sequences[og].append((current_header, ''.join(current_sequence)))
    
    # Write each OG to its own file
    for og, sequences in og_sequences.items():
        output_file = os.path.join(output_directory, f"{og}.fa")
        with open(output_file, 'w') as fasta_file:
            for header, sequence in sequences:
                fasta_file.write(f"{header}\n{sequence}\n")
        print(f"Created {og}.fa with {len(sequences)} sequences")
    
    print(f"\nTotal orthogroups: {len(og_sequences)}")

# ============================================
# Run on Puhti
# ============================================

input_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/bEVs_all_sequences.fa'
output_directory = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/OLD_bEVs_orthogroups_fasta'

split_fasta_by_og(input_fasta, output_directory)