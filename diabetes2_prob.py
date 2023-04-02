import pandas as pd
import numpy as np
import subprocess

# Read in SNP data from CSV file
snp_data = pd.read_csv("snp_data.csv")

# Convert SNP data to PLINK format
plink_data = snp_data.copy()
plink_data.iloc[:,2:] = plink_data.iloc[:,2:].replace({0: "A A", 1: "A B", 2: "B B"})
plink_data.to_csv("plink_data.ped", sep="\t", index=False, header=False)

# Create MAP file
chromosome = "1" # chromosome number
rs_ids = plink_data.iloc[:,0].values # SNP IDs
positions = plink_data.iloc[:,1].values # SNP positions
map_data = np.column_stack((chromosome, rs_ids, 0, positions))
pd.DataFrame(map_data, columns=["CHR", "SNP", "CM", "POS"]).to_csv("plink_data.map", sep="\t", index=False, header=False)

# Convert PLINK data to binary format
plink_cmd = "plink --file plink_data --make-bed --out plink_data"
subprocess.call(plink_cmd, shell=True)


#  use PLINK to perform a logistic regression analysis - diabetes 2 probability?

import subprocess

# Define PLINK command as a string
plink_cmd = "plink --file my_data --pheno diabetes_pheno.txt --logistic --out plink_results"

# Run PLINK command
subprocess.call(plink_cmd, shell=True)

# Read in PLINK results
plink_results = pd.read_csv("plink_results.assoc.logistic", sep="\t")

# Extract relevant statistics
beta = plink_results["BETA"].values[0]
se = plink_results["SE"].values[0]
p_value = plink_results["P"].values[0]

# Calculate odds ratio and probability of developing diabetes
odds_ratio = np.exp(beta)
prob = odds_ratio / (1 + odds_ratio)
