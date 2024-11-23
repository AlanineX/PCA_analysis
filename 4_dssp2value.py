import os
import math
from natsort import natsorted
from multiprocessing import Pool

def parse_dssp_fields(dssp_file_path, extract_fields):
    data = {field: [] for field in extract_fields}
    residue_identifiers = []
    
    with open(dssp_file_path, 'r') as file:
        read_data = False
        for line in file:
            if line.strip().startswith("#  RESIDUE AA STRUCTURE"):
                read_data = True
                continue
            if read_data and len(line) >= 120:
                # Extract residue index and chain ID
                res_index_str = line[5:10].strip()
                chain_id_str = line[11:12].strip()
                res_index = int(res_index_str) if res_index_str else None
                chain_id = chain_id_str if chain_id_str else None
                
                if res_index is not None and chain_id is not None:
                    residue_identifiers.append(f"{res_index}_{chain_id}")

                if "tco" in extract_fields:
                    tco_str = line[85:91].strip()
                    data["tco"].append(float(tco_str) if tco_str else None)

                if "kappa" in extract_fields:
                    kappa_str = line[91:97].strip()
                    data["kappa"].append(float(kappa_str) if kappa_str else None)

                if "alpha" in extract_fields:
                    alpha_str = line[97:103].strip()
                    data["alpha"].append(float(alpha_str) if alpha_str else None)

                if "phi" in extract_fields:
                    phi_str = line[103:109].strip()
                    data["phi"].append(float(phi_str) if phi_str else None)

                if "psi" in extract_fields:
                    psi_str = line[109:115].strip()
                    data["psi"].append(float(psi_str) if psi_str else None)

                if "x" in extract_fields:
                    x_str = line[115:122].strip()
                    data["x"].append(float(x_str) if x_str else None)

                if "y" in extract_fields:
                    y_str = line[122:129].strip()
                    data["y"].append(float(y_str) if y_str else None)

                if "z" in extract_fields:
                    z_str = line[129:136].strip()
                    data["z"].append(float(z_str) if z_str else None)
                    
    return data, residue_identifiers

def transform_data(data, transform):
    for field, values in data.items():
        if transform == "abs":
            data[field] = [abs(value) for value in values if value is not None]
        elif transform == "sin":
            data[field] = [math.sin(math.radians(value)) for value in values if value is not None]
        elif transform == "cos":
            data[field] = [math.cos(math.radians(value)) for value in values if value is not None]
        elif transform == "rad":
            data[field] = [math.radians(value) for value in values if value is not None]
        elif transform == "no":
            data[field] = [value for value in values if value is not None]
    return data

def export_transformed_values(data, residue_identifiers, output_file):
    with open(output_file, 'w') as file:
        header = []
        for res_id in residue_identifiers:
            for field in data.keys():
                header.append(f"{res_id}_{field}")
        file.write(",".join(header) + "\n")

        values_line = []
        for i in range(len(residue_identifiers)):
            for field in data.keys():
                values_line.append(f"{data[field][i]}")
        file.write(",".join(values_line) + "\n")

def process_dssp_file(dssp_file, output_dir, transform="abs", extract_fields=("tco", "kappa", "alpha", "phi", "psi", "x", "y", "z")):
    data, residue_identifiers = parse_dssp_fields(dssp_file, extract_fields)

    transformed_data = transform_data(data, transform)

    base_name = os.path.basename(dssp_file).replace('.dssp', '')
    output_file = os.path.join(output_dir, f"{base_name}_transformed_values.csv")
    export_transformed_values(transformed_data, residue_identifiers, output_file)
    # debug part
    # if "886" in dssp_file:
    #     print(f"Processed {dssp_file} with {len(residue_identifiers)} residues and {len(data)} fields")
    #     print(data)
    #     print(transformed_data)
    #     print(base_name)



def process_all_dssp_files(base_dir, output_dir, num_cores=None, transform="abs", extract_fields=("tco", "kappa", "alpha", "phi", "psi", "x", "y", "z")):
    dssp_files = [os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.endswith('.dssp')]
    os.makedirs(output_dir, exist_ok=True)

    with Pool(processes=num_cores) as pool:
        pool.starmap(process_dssp_file, [(dssp_file, output_dir, transform, extract_fields) for dssp_file in dssp_files])

def aggregate_results_to_csv(output_dir, csv_file, header=None):
    temp_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('_transformed_values.csv')]
    aggregated_data = []

    for temp_file in temp_files:
        with open(temp_file, 'r') as file:
            data = file.readlines()
            base_name = os.path.basename(temp_file).replace('_transformed_values.csv', '')
            data_with_filename = [f"{base_name},{line.strip()}" for line in data[1:]]  
            aggregated_data.append((base_name, data_with_filename))

    aggregated_data = natsorted(aggregated_data, key=lambda x: x[0])

    with open(csv_file, 'w') as file:
        if header:
            file.write(f"filename,{header}\n") 
        for base_name, data_lines in aggregated_data:
            for line in data_lines:
                file.write(line + '\n')

# Change these as needed
base_dir = "./test_monomer/dt1000_fmt_raw"
output_dir = "./test_monomer/dt1000_fmt_raw_value"

num_cores = 2 # Specify the number of cores to use 
transform = "no"  # Options: "abs", "sin", "cos", "rad', "no"
# extract_fields = ("tco", "kappa", "alpha", "phi", "psi", "x", "y", "z")  # Specify fields to extract
extract_fields = ("x", "y", "z")  # for coordinates based PCA
# extract_fields = ("phi", "psi")  # for dihedral angles based PCA

process_all_dssp_files(base_dir, output_dir, num_cores, transform, extract_fields)

sample_dssp_file = os.path.join(base_dir, os.listdir(base_dir)[0])  
_, residue_identifiers = parse_dssp_fields(sample_dssp_file, extract_fields)

header = ",".join([f"{res_id}_{field}" for res_id in residue_identifiers for field in extract_fields])

csv_file = output_dir + "/features_ca.csv" # Change this as needed
aggregate_results_to_csv(output_dir, csv_file, header=header)
