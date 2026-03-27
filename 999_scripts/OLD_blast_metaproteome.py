import pandas as pd
import os
import subprocess

def prepare_epitopes_for_blast(epitope_files, output_dir):
    """
    Create FASTA files of epitopes for BLAST searching.
    """
    
    print("Preparing epitope FASTA files for BLAST...")
    
    fasta_files = {}
    
    for file_name, file_path in epitope_files.items():
        try:
            df = pd.read_csv(file_path, sep='\t')
            
            fasta_file = f"{output_dir}/{file_name}_epitopes.fa"
            
            with open(fasta_file, 'w') as f:
                for _, row in df.iterrows():
                    cluster_id = row['Cluster_ID']
                    sequence = row['Representative_Sequence']
                    
                    f.write(f">{cluster_id}\n")
                    f.write(f"{sequence}\n")
            
            fasta_files[file_name] = fasta_file
            print(f"  Created: {fasta_file}")
            
        except Exception as e:
            print(f"  Error with {file_name}: {e}")
    
    return fasta_files

def blast_epitopes_against_metaproteome(epitope_fasta_files, metaproteome_fasta, output_dir):
    """
    BLAST epitopes against metaproteome to find exact and similar matches.
    """
    
    print("\n" + "="*70)
    print("BLASTING EPITOPES AGAINST METAPROTEOME")
    print("="*70)
    
    # Create BLAST database from metaproteome
    blast_db = f"{output_dir}/metaproteome_blastdb"
    
    print(f"\nCreating BLAST database: {blast_db}")
    cmd = f"makeblastdb -in {metaproteome_fasta} -dbtype prot -out {blast_db}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print("✓ BLAST database created")
    except Exception as e:
        print(f"Error creating BLAST database: {e}")
        return
    
    # BLAST each epitope file
    for file_name, epitope_fasta in epitope_fasta_files.items():
        print(f"\n{'='*70}")
        print(f"BLASTing: {file_name}")
        print(f"{'='*70}")
        
        blast_output = f"{output_dir}/{file_name}_blast_results.txt"
        
        # BLAST parameters optimized for short sequences (epitopes)
        cmd = f"""blastp \
            -query {epitope_fasta} \
            -db {blast_db} \
            -out {blast_output} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
            -task blastp-short \
            -evalue 1000 \
            -word_size 2 \
            -max_target_seqs 1000 \
            -num_threads 4
        """
        
        print(f"Running BLAST (this may take a few minutes)...")
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"✓ BLAST complete: {blast_output}")
        except Exception as e:
            print(f"Error running BLAST: {e}")
    
    print("\n" + "="*70)
    print("BLAST COMPLETE")
    print("="*70)

def parse_blast_results(epitope_files, blast_dir, metaproteome_fasta, output_dir):
    """
    Parse BLAST results and create comprehensive summary.
    Fixed to exclude LITERATURE from genera lists.
    """
    
    print("\n" + "="*70)
    print("PARSING BLAST RESULTS")
    print("="*70)
    
    # Load metaproteome metadata
    print("\nLoading metaproteome metadata...")
    metaproteome_meta = {}
    
    with open(metaproteome_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Parse: >ProteinID|Intensity|Genus|Gram|Name
                parts = line[1:].strip().split('|')
                protein_id = parts[0]
                metaproteome_meta[protein_id] = {
                    'intensity': float(parts[1]) if len(parts) > 1 else 0,
                    'genus': parts[2] if len(parts) > 2 else 'Unknown',
                    'gram': parts[3] if len(parts) > 3 else 'Unknown',
                    'name': parts[4] if len(parts) > 4 else 'Unknown'
                }
    
    print(f"Loaded metadata for {len(metaproteome_meta)} proteins")
    
    # Parse each BLAST result
    for file_name, epitope_file_path in epitope_files.items():
        print(f"\n{'='*70}")
        print(f"Processing: {file_name}")
        print(f"{'='*70}")
        
        # Load epitope metadata
        try:
            epitope_df = pd.read_csv(epitope_file_path, sep='\t')
            epitope_meta = {}
            for _, row in epitope_df.iterrows():
                epitope_meta[row['Cluster_ID']] = row.to_dict()
        except Exception as e:
            print(f"Error loading epitope file: {e}")
            continue
        
        # Load BLAST results
        blast_file = f"{blast_dir}/{file_name}_blast_results.txt"
        
        try:
            blast_df = pd.read_csv(blast_file, sep='\t', header=None, names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq'
            ])
            print(f"Loaded {len(blast_df)} BLAST hits")
        except Exception as e:
            print(f"No BLAST results for {file_name}: {e}")
            continue
        
        # Aggregate by epitope
        results = []
        
        for cluster_id, meta in epitope_meta.items():
            epitope_seq = meta['Representative_Sequence']
            
            # Get all hits for this epitope
            hits = blast_df[blast_df['qseqid'] == cluster_id]
            
            # Categorize hits by identity
            exact_matches = hits[hits['pident'] == 100]
            high_similarity = hits[(hits['pident'] >= 90) & (hits['pident'] < 100)]
            medium_similarity = hits[(hits['pident'] >= 80) & (hits['pident'] < 90)]
            low_similarity = hits[(hits['pident'] >= 70) & (hits['pident'] < 80)]
            
            # Calculate total intensity for each category
            def get_total_intensity(hit_df):
                if len(hit_df) == 0:
                    return 0
                total = 0
                for protein_id in hit_df['sseqid'].unique():
                    total += metaproteome_meta.get(protein_id, {}).get('intensity', 0)
                return total
            
            exact_intensity = get_total_intensity(exact_matches)
            high_sim_intensity = get_total_intensity(high_similarity)
            medium_sim_intensity = get_total_intensity(medium_similarity)
            low_sim_intensity = get_total_intensity(low_similarity)
            
            # Get genera for matches - FIXED to exclude LITERATURE and Unknown
            def get_genera(hit_df):
                if len(hit_df) == 0:
                    return 'None'
                genera = set()
                for protein_id in hit_df['sseqid'].unique():
                    genus = metaproteome_meta.get(protein_id, {}).get('genus', 'Unknown')
                    # FILTER OUT LITERATURE and Unknown
                    if genus not in ['LITERATURE', 'Unknown', '']:
                        genera.add(genus)
                return '; '.join(sorted(genera)) if genera else 'None'
            
            # Get example protein IDs
            def get_protein_examples(hit_df, max_examples=5):
                if len(hit_df) == 0:
                    return 'None'
                proteins = hit_df['sseqid'].unique()[:max_examples]
                return '; '.join(proteins)
            
            results.append({
                'Cluster_ID': cluster_id,
                'Rank': meta.get('Rank', ''),
                'Epitope_Sequence': epitope_seq,
                'Sequence_Length': len(epitope_seq),
                'Priority_Score': meta.get('Priority_Score', ''),
                'Methods_Used': meta.get('Methods_Used', ''),
                'IEDB_Links': meta.get('IEDB_Links', ''),
                'Predicted_Species': meta.get('Species_List', ''),
                'Predicted_Num_Species': meta.get('Num_Species', ''),
                'Suitable_For_Fishing': meta.get('Suitable_For_Fishing', ''),
                
                # Exact matches (100% identity)
                'Exact_Match_Count': len(exact_matches),
                'Exact_Match_Proteins': len(exact_matches['sseqid'].unique()) if len(exact_matches) > 0 else 0,
                'Exact_Match_Intensity': exact_intensity,
                'Exact_Match_Genera': get_genera(exact_matches),
                'Exact_Match_Example_Proteins': get_protein_examples(exact_matches),
                
                # High similarity (90-99%)
                'High_Similarity_Count': len(high_similarity),
                'High_Similarity_Proteins': len(high_similarity['sseqid'].unique()) if len(high_similarity) > 0 else 0,
                'High_Similarity_Intensity': high_sim_intensity,
                'High_Similarity_Genera': get_genera(high_similarity),
                'High_Similarity_Example_Proteins': get_protein_examples(high_similarity),
                
                # Medium similarity (80-89%)
                'Medium_Similarity_Count': len(medium_similarity),
                'Medium_Similarity_Proteins': len(medium_similarity['sseqid'].unique()) if len(medium_similarity) > 0 else 0,
                'Medium_Similarity_Intensity': medium_sim_intensity,
                'Medium_Similarity_Genera': get_genera(medium_similarity),
                
                # Low similarity (70-79%)
                'Low_Similarity_Count': len(low_similarity),
                'Low_Similarity_Proteins': len(low_similarity['sseqid'].unique()) if len(low_similarity) > 0 else 0,
                'Low_Similarity_Intensity': low_sim_intensity,
                
                # Total
                'Total_Hits': len(hits),
                'Total_Unique_Proteins': len(hits['sseqid'].unique()) if len(hits) > 0 else 0,
                'Total_Intensity': exact_intensity + high_sim_intensity + medium_sim_intensity + low_sim_intensity,
                'Any_Match_Found': 'Yes' if len(hits) > 0 else 'No'
            })
        
        # Create DataFrame
        results_df = pd.DataFrame(results)
        
        # Sort by exact matches first, then high similarity, then total intensity
        results_df = results_df.sort_values(
            ['Exact_Match_Count', 'High_Similarity_Count', 'Total_Intensity'],
            ascending=[False, False, False]
        )
        
        # Save
        output_file = f"{output_dir}/{file_name}_blast_summary.txt"
        results_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"✓ Saved: {output_file}")
        
        # Statistics
        with_exact = len(results_df[results_df['Exact_Match_Count'] > 0])
        with_high_sim = len(results_df[results_df['High_Similarity_Count'] > 0])
        with_medium_sim = len(results_df[results_df['Medium_Similarity_Count'] > 0])
        with_any = len(results_df[results_df['Any_Match_Found'] == 'Yes'])
        total = len(results_df)
        
        print(f"\nResults:")
        print(f"  Exact matches (100%): {with_exact}/{total} ({100*with_exact/total:.1f}%)")
        print(f"  High similarity (90-99%): {with_high_sim}/{total} ({100*with_high_sim/total:.1f}%)")
        print(f"  Medium similarity (80-89%): {with_medium_sim}/{total} ({100*with_medium_sim/total:.1f}%)")
        print(f"  Any match (≥70%): {with_any}/{total} ({100*with_any/total:.1f}%)")
        
        if with_exact > 0:
            print(f"\nTop 5 exact matches by intensity:")
            top_exact = results_df[results_df['Exact_Match_Count'] > 0].head(5)
            for _, row in top_exact.iterrows():
                print(f"  {row['Epitope_Sequence'][:20]:20} - {row['Exact_Match_Proteins']} proteins, intensity={row['Exact_Match_Intensity']:.2e}")

print("\n" + "="*70)
print("BLAST ANALYSIS COMPLETE")
print("="*70)
print("\nOutput files:")
print("  *_blast_summary.txt - Aggregated results by epitope")
print("  *_blast_results.txt - Raw BLAST output")
print("="*70)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

cluster_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
metaproteome_fasta = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/001_data/metaproteome_sequences.fa'
output_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'

# Epitope files to analyze
epitope_files = {
    'TOP50_FINAL': f"{cluster_dir}/TOP50_FINAL_with_metadata.txt",
    'COMBINED_TOP_TIERS_FINAL': f"{cluster_dir}/COMBINED_TOP_TIERS_FINAL_with_metadata.txt",
    'COMBINED_FISHABLE': f"{cluster_dir}/COMBINED_TOP_TIERS_FISHABLE_epitopes.txt"
}

print("="*70)
print("EPITOPE-METAPROTEOME BLAST ANALYSIS")
print("="*70)
print(f"\nMetaproteome: {metaproteome_fasta}")
print(f"Output directory: {output_dir}")
print(f"Epitope files to process: {len(epitope_files)}")
print("="*70)

# Step 1: Prepare epitope FASTA files
fasta_files = prepare_epitopes_for_blast(epitope_files, output_dir)

# Step 2: Run BLAST
blast_epitopes_against_metaproteome(fasta_files, metaproteome_fasta, output_dir)

# Step 3: Parse and summarize results
parse_blast_results(epitope_files, output_dir, metaproteome_fasta, output_dir)