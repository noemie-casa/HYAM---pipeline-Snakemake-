#HYAM_stage 
#Missing data for SNP

################################################## 

#Packages to load 
import pandas as pd 
import numpy as np 
import sys 
import matplotlib.pyplot as plt 

# Implémentation des arguments 
missing_file = sys.argv[1] 
output_file = sys.argv[2]

df = pd.read_csv(missing_file, sep="\t")
missing_data = df["F_MISS"]

# plot histogram 
plt.figure() 
plt.hist(missing_data, bins=20) 
plt.xlabel("Proportion of missing genotypes") 
plt.ylabel("Number of SNPs") 
plt.title("SNP missingness distribution") 
plt.savefig(output_file) 
