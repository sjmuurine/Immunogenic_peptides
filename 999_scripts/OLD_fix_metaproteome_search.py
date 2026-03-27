import pandas as pd
import os

def fix_metaproteome_search_files(search_dir):
    """
    Remove 'LITERATURE' from Predicted_Species_List in metaproteome search files.
    """
    
    print("="*70)
    print("FIXING METAPROTEOME SEARCH FILES")
    print("="*70)
    
    # Find all metaproteome search result files
    files = [f for f in os.listdir(search_dir) if f.endswith('_metaproteome_search.txt')]
    
    for filename in files:
        filepath = os.path.join(search_dir, filename)
        
        print(f"\nProcessing: {filename}")
        
        try:
            df = pd.read_csv(filepath, sep='\t')
            
            # Fix Predicted_Species_List column
            def clean_species_list(species_str):
                if pd.isna(species_str) or species_str == 'None':
                    return species_str
                
                # Split, filter out LITERATURE, rejoin
                species_list = [s.strip() for s in str(species_str).split(';')]
                cleaned = [s for s in species_list if s not in ['LITERATURE', 'Unknown', '']]
                
                return '; '.join(sorted(set(cleaned))) if cleaned else 'None'
            
            # Apply fix
            if 'Predicted_Species_List' in df.columns:
                original = df['Predicted_Species_List'].copy()
                df['Predicted_Species_List'] = df['Predicted_Species_List'].apply(clean_species_list)
                
                # Update Predicted_Num_Species
                if 'Predicted_Num_Species' in df.columns:
                    df['Predicted_Num_Species'] = df['Predicted_Species_List'].apply(
                        lambda x: len([s for s in str(x).split(';') if s.strip() and s.strip() != 'None']) if pd.notna(x) and x != 'None' else 0
                    )
                
                changed = sum(original != df['Predicted_Species_List'])
                print(f"  Fixed {changed} rows in Predicted_Species_List")
            
            # Save
            df.to_csv(filepath, sep='\t', index=False)
            print(f"  ✓ Saved")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "="*70)
    print("FIXING COMPLETE")
    print("="*70)

# Run
search_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'
fix_metaproteome_search_files(search_dir)