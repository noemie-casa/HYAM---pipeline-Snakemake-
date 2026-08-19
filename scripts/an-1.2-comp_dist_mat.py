#HYAM_stage 
#Distance matrix for all chromosomes 

################################################## 

#Packages to load 
import pandas as pd 
import numpy as np 
import sys 

# Implémentation des arguments 
filtered_file = sys.argv[1] 
output_file = sys.argv[2] 

df = pd.read_csv(filtered_file, sep="\t")
chromosomes = df.iloc[:,0].values
positions = df.iloc[:,1].values 

# distances 
distances = []
for chrom in np.unique(chromosomes):
    idx = np.where(chromosomes == chrom)[0]
    pos = positions[idx]
    dist = np.diff(pos)
    for i in range (len(dist)):
        distances.append([chrom, pos[i], pos[i+1], dist[i]])

# Distance matrix construction
pd.DataFrame(distances).to_csv(output_file, sep="\t", index=False)
