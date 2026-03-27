import pandas as pd

# Quick test of parsing
test_header = "WP_239447903.1 V_magna_OG0000080_gram+_Cytoplasmic Membrane"

species_map = {
    'A': 'A_muciniphila',
    'B': 'B_fragilis',
    'E': 'E_cloacae',
    'F': 'F_prausnitzii',
    'L': 'L_reuteri',
    'P': 'P_russellii',
    'S': 'S_salivarius',
    'T': 'T_sanguinis',
    'V': 'V_magna'
}

# Split by space
parts = test_header.split(' ', 1)
print(f"After space split: {parts}")

protein_id = parts[0]
metadata = parts[1]
print(f"Protein ID: {protein_id}")
print(f"Metadata: {metadata}")

# Split metadata by underscore
meta_parts = metadata.split('_')
print(f"After underscore split: {meta_parts}")
print(f"Length: {len(meta_parts)}")

if len(meta_parts) >= 4:
    species_abbrev = meta_parts[0]
    orthogroup = meta_parts[1]
    conservation = meta_parts[2]
    localization = ' '.join(meta_parts[3:])
    
    print(f"\nExtracted:")
    print(f"  Species abbrev: '{species_abbrev}'")
    print(f"  Orthogroup: '{orthogroup}'")
    print(f"  Conservation: '{conservation}'")
    print(f"  Localization: '{localization}'")
    
    species_full = species_map.get(species_abbrev, species_abbrev)
    print(f"  Species full: '{species_full}'")