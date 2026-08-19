#HYAM_stage 
#Distance distribution plot for all chromosomes

################################################## 

#Packages to load 
import pandas as pd 
import numpy as np 
import sys 
import matplotlib.pyplot as plt 

# Implémentation des arguments 
dist_file = sys.argv[1] 
output_file = sys.argv[2] 
method = sys.argv[3] 

df = pd.read_csv(dist_file, sep="\t") 
distances = df.iloc[1:, 3] 

# plot histogram 
plt.figure() 
plt.hist(np.log10(distances), bins=20) 
plt.xlabel("log10(Distance (b))") 
plt.ylabel("Number of loci") 
plt.title("Distance distribution of loci for the " + method + " method") 
plt.savefig(output_file)


