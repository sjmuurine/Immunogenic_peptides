import pandas as pd
import os

def fix_species_list_literature(metadata_dir):
    """
    Remove 'LITERATURE' from Species_List column in metadata files.
    """
    
    print("="*70)
    print("FIXING SPECIES LISTS IN METADATA FILES")
    print("="*70)
    
    # Files to fix
    files_to_fix = [
        'TOP50_FINAL_with_metadata.txt',
        'TIER1_FINAL_with_metadata.txt',
        'TIER2_FINAL_with_metadata.txt',
        'COMBINED_TOP_TIERS_FINAL_with_metadata.txt',
        'TOP50_FISHABLE_epitopes.txt',
        'TIER1_FISHABLE_epitopes.txt',
        'TIER2_FISHABLE_epitopes.txt',
        'COMBINED_TOP_TIERS_FISHABLE_epitopes.txt'
    ]
    
    for filename in files_to_fix:
        filepath = os.path.join(metadata_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"\nSkipping {filename} (not found)")
            continue
        
        print(f"\nProcessing: {filename}")
        
        try:
            df = pd.read_csv(filepath, sep='\t')
            
            # Fix Species_List column
            def clean_species_list(species_str):
                if pd.isna(species_str) or species_str == 'None':
                    return species_str
                
                # Split by semicolon, filter out LITERATURE, rejoin
                species_list = [s.strip() for s in str(species_str).split(';')]
                cleaned = [s for s in species_list if s not in ['LITERATURE', 'Unknown', '']]
                
                # Update Num_Species if column exists
                return '; '.join(sorted(set(cleaned))) if cleaned else 'None'
            
            # Apply fix to Species_List
            if 'Species_List' in df.columns:
                original = df['Species_List'].copy()
                df['Species_List'] = df['Species_List'].apply(clean_species_list)
                
                # Update Num_Species to match
                if 'Num_Species' in df.columns:
                    df['Num_Species'] = df['Species_List'].apply(
                        lambda x: len([s for s in str(x).split(';') if s.strip() and s.strip() != 'None']) if pd.notna(x) and x != 'None' else 0
                    )
                
                # Count changes
                changed = sum(original != df['Species_List'])
                print(f"  Fixed {changed} rows")
            else:
                print(f"  No Species_List column found")
            
            # Also fix Predicted_Species_List if it exists (in metaproteome search files)
            if 'Predicted_Species_List' in df.columns:
                df['Predicted_Species_List'] = df['Predicted_Species_List'].apply(clean_species_list)
                
                if 'Predicted_Num_Species' in df.columns:
                    df['Predicted_Num_Species'] = df['Predicted_Species_List'].apply(
                        lambda x: len([s for s in str(x).split(';') if s.strip() and s.strip() != 'None']) if pd.notna(x) and x != 'None' else 0
                    )
            
            # Save
            df.to_csv(filepath, sep='\t', index=False)
            print(f"  ✓ Saved")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "="*70)
    print("FIXING COMPLETE")
    print("="*70)

# Run on metadata files
metadata_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/004_cluster'
fix_species_list_literature(metadata_dir)

# Also run on metaproteome search results
metaproteome_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'
fix_species_list_literature(metaproteome_dir)