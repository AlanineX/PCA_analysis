import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

def run_dssp(input_file, output_file, log_file_path):
    # command = f"dssp -i {input_file} -o {output_file}"
    command = f"mkdssp {input_file} {output_file}"
    try:
        subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Processed {input_file}")
    except subprocess.CalledProcessError as e:
        error_message = f"Error processing {input_file}: {e.stderr.decode()}\n"
        print(error_message)
        with open(log_file_path, 'a') as log_file:
            log_file.write(error_message)

def main():
    input_dir = "./test_monomer/dt1000_fmt"
    output_dir = "./test_monomer/dt1000_fmt_raw"
    log_file_path = "./dssp_errors.log"
    
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = [f for f in os.listdir(input_dir) if f.endswith('.pdb')]
    
    # do not use all cores by os.cpu_count()
    with ThreadPoolExecutor(max_workers=4) as executor:
        for pdb_file in pdb_files:
            input_file = os.path.join(input_dir, pdb_file)
            output_file = os.path.join(output_dir, pdb_file.replace('.pdb', '.dssp'))
            executor.submit(run_dssp, input_file, output_file, log_file_path)

if __name__ == "__main__":
    main()
