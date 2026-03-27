import csv
import os

def separate_fasta_by_orthogroup(input_file, output_directory):
    """
    Creates separate FASTA files for each orthogroup.
    
    Organizes sequences into:
    output_directory/
        OG0000089.fa
        OG0000090.fa
        ...
    
    Skips sequences with orthogroup = "unique"
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    
    # Dictionary to store sequences by OG
    og_sequences = {}
    
    with open(input_file, 'r', encoding='utf-8') as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter='\t')
        
        for row in reader:
            # Extract components
            orthogroup = row['Orthogroup'].strip()
            conserved_in = row['Conserved in'].strip()
            fasta_headers = row['Fasta headers'].strip()
            localization = row['Localization'].strip()
            sequence = row['Sequence (primary major protein)'].strip()
            
            # Skip if no sequence or if "unique"
            if not sequence or orthogroup.lower() == "unique" or not orthogroup:
                continue
            
            # Handle missing values
            if not conserved_in:
                conserved_in = "NA"
            if not localization:
                localization = "Unknown"
            
            # Create header
            header = f"{orthogroup}_{conserved_in}_{fasta_headers}_{localization}"
            
            # Add to dictionary
            if orthogroup not in og_sequences:
                og_sequences[orthogroup] = []
            og_sequences[orthogroup].append((header, sequence))
    
    # Write each OG to its own file
    for og, sequences in og_sequences.items():
        output_file = os.path.join(output_directory, f"{og}.fa")
        with open(output_file, 'w', encoding='utf-8') as fasta_file:
            for header, sequence in sequences:
                fasta_file.write(f">{header}\n{sequence}\n")
        print(f"Created {output_file} with {len(sequences)} sequences")
    
    print(f"\nTotal orthogroups: {len(og_sequences)}")

# ============================================
# Run the script in your local computer using TSV file from the excel having results
# ============================================

input_file = '/path/to/your/big_table.tsv'
output_directory = '/path/to/output_orthogroups/'

separate_fasta_by_orthogroup(input_file, output_directory)