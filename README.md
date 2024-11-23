# PCA_analysis mPro_example_1122

Install the necessary Python packages:
```bash
pip install numpy matplotlib scikit-learn biopython mdanalysis pandas
```

Additional Software
DSSP: Required for secondary structure processing. Install DSSP using your system's package manager (e.g., ```bash sudo apt install dssp ```)

## File Descriptions
- **`1_pdb_split.py`**: Splits a single trajectory PDB file into individual files.
- **`2_pdb_DSSP_format.py`**: Converts PDB files into a format recognized by DSSP.
- **`3_pdb2dssp.py`**: Processes PDB files and extracts DSSP raw data.
- **`4_dssp2value.py`**: Extracts and transforms DSSP information.
- **`5_kmeans_PCA.py`**: Performs PCA and K-means clustering.
