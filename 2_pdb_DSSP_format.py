from Bio import PDB
import os
import shutil

def reformat_pdb(input_file, output_file, chain_id='A'):
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_file)

    # # Set chain identifier for all chains in the structure
    # for model in structure:
    #     for chain in model:
    #         chain.id = chain_id

    # Dictionary to map non-standard residue names to standard ones
    residue_name_mapping = {
        'HSP': 'HIS',
        'GLX': 'GLU',
        'ASX': 'ASP',
        'D6M': 'LYS',
        'DMH': 'LYS'
    }

    # Replace specific non-standard residues with their standard counterparts
    for model in structure:
        for chain in model:
            # format_c_term_oxygen(chain)  # Call the function to format C-terminal oxygen, not work when saving file
            for residue in chain:
                if residue.resname in residue_name_mapping:
                    residue.resname = residue_name_mapping[residue.resname]
                    print(f"renaming{residue.resname}")

    
    # Save the reformatted structure to the output file
    io = PDB.PDBIO()
    io.set_structure(structure)
    
    with open(output_file, 'w') as out:
        out.write("HEADER    EXAMPLE STRUCTURE\n") # Write a HEADER line
        # Write a generic CRYST1 line
        out.write("CRYST1   1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n")
        io.save(out)

def replace_ot_in_file(input_file, output_file):
    with open(input_file, 'r') as f:
        content = f.read()

    # Replace OT1 and OT2 with O and OXT respectively
    content = content.replace(' OT1 ', ' O   ')
    content = content.replace(' OT2 ', ' OXT ')

    with open(output_file, 'w') as f:
        f.write(content)
        # f.write('   4.95219   4.95219   4.95219')

base_dir = './test_monomer/dt1000'
formatted_dir = './test_monomer/dt1000_fmt'

# Create the new directory if it doesn't exist
if not os.path.exists(formatted_dir):
    os.makedirs(formatted_dir)


for i in range(1, 1002):
    input_file_name = f'{i}.pdb'
    'revise here if needed'
    # input_file_name = '1100_1200_sys.pdb'
    input_pdb_path = os.path.join(base_dir, input_file_name)
    temp_output_path = os.path.join(base_dir, input_file_name.replace('.pdb', '_temp.pdb'))
    final_output_path = os.path.join(formatted_dir, input_file_name.replace('.pdb', '_fmt.pdb'))

    if os.path.isfile(input_pdb_path):
        reformat_pdb(input_pdb_path, temp_output_path)
        replace_ot_in_file(temp_output_path, final_output_path)

        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)

        # shutil.move(final_output_path, os.path.join(formatted_dir, os.path.basename(final_output_path)))
        print(f'Processing {input_file_name} ...')
    else:
        print(f'File {input_file_name} does not exist. Skipping.')

print('All done.')
