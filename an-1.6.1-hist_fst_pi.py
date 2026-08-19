# HYAM_stage 
# Fst and Pi distribuition between two species 

################################################## 

#Packages to load 
import pandas as pd 
import numpy as np 
import sys 
import matplotlib.pyplot as plt 

# Implémentation des arguments 
Fst_file = sys.argv[1] 
Pi_fileA = sys.argv[2] 
Pi_fileB = sys.argv[3]
output_file = sys.argv[4]
speciesA = sys.argv[5]
speciesB = sys.argv[6]

Fst_df = pd.read_csv(Fst_file, sep="\t") 
PiA_df = pd.read_csv(Pi_fileA, sep="\t")
PiB_df = pd.read_csv(Pi_fileB, sep="\t")

Fst = Fst_df.iloc[:,2].dropna()
PiA = PiA_df.iloc[:,2].dropna()
PiB = PiB_df.iloc[:,2].dropna()

# Moyenne et Seuil
mean_Fst = np.mean(Fst)
q95_Fst = np.quantile(Fst, 0.95)

mean_PiA = np.mean(PiA)
q95_PiA = np.quantile(PiA, 0.95)

mean_PiB = np.mean(PiB)
q95_PiB = np.quantile(PiB, 0.95)

# plot histograms 
fig, axes = plt.subplots(3, 1, figsize=(8, 12))

# Boxes Pi 
all_pi = np.concatenate([PiA, PiB])
bins = np.linspace(all_pi.min(), all_pi.max(), 21)

# Fst
axes[0].hist(Fst, bins=20)
axes[0].axvline(mean_Fst, linestyle="--", color="green", label=f"Mean = {mean_Fst:.4f}")
axes[0].axvline(q95_Fst, linestyle="--", color="red", label=f"95% quantile = {q95_Fst:.4f}")
axes[0].set_title(f"Fst distribution between {speciesA} and {speciesB}")
axes[0].set_xlabel("Fst")
axes[0].set_ylabel("Number of windows")
axes[0].legend()

# PiA
axes[1].hist(PiA, bins=bins)
axes[1].axvline(mean_PiA, linestyle="--", color="green", label=f"Mean = {mean_PiA:.4f}")
axes[1].axvline(q95_PiA, linestyle="--", color="red", label=f"95% quantile = {q95_PiA:.4f}")
axes[1].set_xlim(all_pi.min(), all_pi.max())
axes[1].set_title(f"Pi distribution - {speciesA}")
axes[1].set_xlabel("Pi")
axes[1].set_ylabel("Number of windows")
axes[1].legend()

# PiB
axes[2].hist(PiB, bins=bins)
axes[2].axvline(mean_PiB, linestyle="--", color="green", label=f"Mean = {mean_PiB:.4f}")
axes[2].axvline(q95_PiB, linestyle="--", color="red", label=f"95% quantile = {q95_PiB:.4f}")
axes[2].set_xlim(all_pi.min(), all_pi.max())
axes[2].set_title(f"Pi distribution - {speciesB}")
axes[2].set_xlabel("Pi")
axes[2].set_ylabel("Number of windows")
axes[2].legend()

plt.tight_layout()
plt.savefig(output_file)
plt.close()


