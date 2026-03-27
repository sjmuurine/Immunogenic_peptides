import os

def add_blast_headers(blast_dir):
    """
    Add headers to raw BLAST output files.
    """
    
    print("Adding headers to BLAST files...")
    
    # BLAST -outfmt 6 column headers
    headers = [
        'Query_ID',
        'Subject_ID', 
        'Percent_Identity',
        'Alignment_Length',
        'Mismatches',
        'Gap_Opens',
        'Query_Start',
        'Query_End',
        'Subject_Start',
        'Subject_End',
        'E_value',
        'Bit_Score',
        'Query_Sequence',
        'Subject_Sequence'
    ]
    
    header_line = '\t'.join(headers) + '\n'
    
    # Find BLAST result files (not summaries)
    files = [f for f in os.listdir(blast_dir) if f.endswith('_blast_results.txt')]
    
    for filename in files:
        filepath = os.path.join(blast_dir, filename)
        
        print(f"Processing: {filename}")
        
        # Read existing content
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Write with header
        with open(filepath, 'w') as f:
            f.write(header_line)
            f.write(content)
        
        print(f"  ✓ Header added")
    
    print("\nComplete!")

# Run
blast_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'
add_blast_headers(blast_dir)