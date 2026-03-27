import pandas as pd
import os

output_dir = "/scratch/project_2009813/JOHANNA/II_potential_epitopes/003_result_comparisons"
temp_dir = f"{output_dir}/temp"

# Create temp directory
os.makedirs(temp_dir, exist_ok=True)

# Collect all unique sequences
all_sequences = set()
sequence_sources = {}  # Track which source each sequence came from

# Load HIGH_CONFIDENCE overlapping epitopes
print("Loading HIGH_CONFIDENCE epitopes...")
high_conf = pd.read_csv(f"{output_dir}/HIGH_CONFIDENCE_overlapping_epitopes.txt", sep='\t')
for seq in high_conf['Sequence']:
    all_sequences.add(seq)
    if seq not in sequence_sources:
        sequence_sources[seq] = []
    sequence_sources[seq].append('HIGH_CONFIDENCE')

# Load IEDB_only
print("Loading IEDB_only epitopes...")
iedb_only = pd.read_csv(f"{output_dir}/IEDB_only_epitopes.txt", sep='\t')
for seq in iedb_only['Sequence']:
    all_sequences.add(seq)
    if seq not in sequence_sources:
        sequence_sources[seq] = []
    sequence_sources[seq].append('IEDB_only')

# Load Kmer_EXACT_only
print("Loading Kmer_EXACT_only epitopes...")
kmer_exact = pd.read_csv(f"{output_dir}/Kmer_EXACT_only_epitopes.txt", sep='\t')
for seq in kmer_exact['Sequence']:
    all_sequences.add(seq)
    if seq not in sequence_sources:
        sequence_sources[seq] = []
    sequence_sources[seq].append('Kmer_EXACT_only')

# Load Kmer_FUZZY_only
print("Loading Kmer_FUZZY_only epitopes...")
kmer_fuzzy = pd.read_csv(f"{output_dir}/Kmer_FUZZY_only_epitopes.txt", sep='\t')
for seq in kmer_fuzzy['Sequence']:
    all_sequences.add(seq)
    if seq not in sequence_sources:
        sequence_sources[seq] = []
    sequence_sources[seq].append('Kmer_FUZZY_only')

print(f"\nTotal unique epitope sequences: {len(all_sequences)}")

# Write to FASTA
fasta_file = f"{temp_dir}/all_epitopes_combined.fasta"
with open(fasta_file, 'w') as f:
    for i, seq in enumerate(sorted(all_sequences), 1):
        sources = ','.join(sequence_sources[seq])
        f.write(f">epitope_{i}|{sources}|len_{len(seq)}\n{seq}\n")

print(f"✓ FASTA created: {fasta_file}")

# Save sequence-to-source mapping
import csv
with open(f"{temp_dir}/epitope_source_mapping.txt", 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Epitope_ID', 'Sequence', 'Length', 'Source_Methods'])
    for i, seq in enumerate(sorted(all_sequences), 1):
        sources = ','.join(sequence_sources[seq])
        writer.writerow([f"epitope_{i}", seq, len(seq), sources])

print("✓ Source mapping saved")
print(f"✓ Total epitopes: {len(all_sequences)}")