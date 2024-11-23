#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 15:56:39 2023

@author: alan
"""

from Bio import PDB
import os

def split_pdb(input_pdb, output_directory="."):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)
    name_number = 0
    
    for model in structure:
        # skip the first one or no
        if model.id == 0:
            continue
        
        output_file = f"{output_directory}/{name_number + model.id}.pdb"
        io = PDB.PDBIO()
        io.set_structure(model)
        io.save(output_file)
        print(f"Saved model {name_number + model.id} to {output_file}")

if __name__ == "__main__":
    input_file = "./test_monomer/7jpy_n17.pdb"
    output_dir = "./test_monomer/dt1000"  
    os.makedirs(output_dir, exist_ok=True)
    split_pdb(input_file, output_dir)
