#!/usr/bin/env python3
import os
import re

# List of missing files found in the warning messages
missing_files = [
    # Input files
    'input_files/fort13.rst',
    'input_files/fort19.rst',
    'input_files/fort20.rst',
    'input_files/fort22.rst',
    'input_files/fort23.rst',
    'input_files/fort24.rst',
    'input_files/meteorological_inputs.rst',
    'input_files/river_inputs.rst',
    'input_files/scalar_inputs.rst',
    'input_files/bathymetry_inputs.rst',
    'input_files/hotstart.rst',
    
    # Output files
    'output_files/velocity_outputs.rst',
    'output_files/harmonic_outputs.rst',
    'output_files/meteorological_outputs.rst',
    'output_files/max_min_outputs.rst',
    'output_files/3d_outputs.rst',
    'output_files/inundation_outputs.rst',
    'output_files/additional_outputs.rst',
    
    # Parameters
    'parameters/nodal_attributes.rst',
    'parameters/boundary_conditions.rst',
    
    # Version history
    'version_history/previous_versions.rst',
]

def create_stub_file(file_path):
    # Get the base filename without extension
    base_name = os.path.basename(file_path).replace('.rst', '')
    
    # Convert snake_case to Title Case
    title = ' '.join(word.capitalize() for word in base_name.split('_'))
    
    # Special case for numbered files like fort13, fort19, etc.
    if base_name.startswith('fort'):
        number = base_name[4:]
        title = f"Fort.{number} File"
    elif base_name == '3d_outputs':
        title = "3D Output Files"
    
    # Create content for the stub file
    content = f"{title}\n{'=' * len(title)}\n\n"
    content += ".. note::\n"
    content += "   This page is under development. Content will be added soon.\n"
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Write the file
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"Created: {file_path}")

def main():
    docs_dir = os.path.dirname(os.path.abspath(__file__))
    
    for file_rel_path in missing_files:
        file_path = os.path.join(docs_dir, file_rel_path)
        create_stub_file(file_path)
    
    print(f"Created {len(missing_files)} stub files.")

if __name__ == "__main__":
    main() 