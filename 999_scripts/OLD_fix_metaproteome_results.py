import pandas as pd
import os

def fix_metaproteome_search_results(results_dir):
    """
    Fix LITERATURE appearing in species/genera lists.
    """
    
    print("="*70)
    print("FIXING METAPROTEOME SEARCH RESULTS")
    print("="*70)
    
    # Find all metaproteome search result files
    files = [f for f in os.listdir(results_dir) if f.endswith('_metaproteome_search.txt')]
    
    for filename in files:
        filepath = os.path.join(results_dir, filename)
        
        print(f"\nProcessing: {filename}")
        
        try:
            df = pd.read_csv(filepath, sep='\t')
            
            # Fix Metaproteome_Genera column
            def clean_genera(genera_str):
                if pd.isna(genera_str) or genera_str == 'None':
                    return genera_str
                
                # Split, filter out LITERATURE and Unknown, rejoin
                genera_list = [g.strip() for g in str(genera_str).split(';')]
                cleaned = [g for g in genera_list if g not in ['LITERATURE', 'Unknown', '']]
                
                return '; '.join(sorted(cleaned)) if cleaned else 'None'
            
            # Apply fix
            if 'Metaproteome_Genera' in df.columns:
                df['Metaproteome_Genera'] = df['Metaproteome_Genera'].apply(clean_genera)
            
            # Save
            df.to_csv(filepath, sep='\t', index=False)
            print(f"  ✓ Fixed")
            
        except Exception as e:
            print(f"  Error: {e}")
    
    print("\n" + "="*70)
    print("FIXING COMPLETE")
    print("="*70)

# Run
results_dir = '/scratch/project_2009813/JOHANNA/II_potential_epitopes/005_metaproteome_validation'
fix_metaproteome_search_results(results_dir)