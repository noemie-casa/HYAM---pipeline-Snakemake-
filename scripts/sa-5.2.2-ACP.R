#HYAM_stage
#ACP et DAPC

##################################################

#Packages to load
library(ggplot2)
library(adegenet)
library(ggrepel)
library(vcfR)

# Files 
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
outdir <- args[2]
nf <- as.numeric(args[3])
method <- args[4]
species_file <- args[5]

# Import vcf
vcf <- read.vcfR(vcf_file)
genlight <- vcfR2genlight(vcf)
species <- read.table(species_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# PCA
pca_res <- glPca(genlight, nf = nf, loadings = TRUE)
percent_pca <- round(pca_res$eig / sum(pca_res$eig) * 100, 2)
names(percent_pca) <- paste0("PC", seq_along(percent_pca))

# Individuals coordinates
scores <- as.data.frame(pca_res$scores)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
scores$INDIVIDUALS <- rownames(scores)

scores <- merge(scores, species, by.x = "INDIVIDUALS", all.x = TRUE)

# PCA plots
pairs <- list()
for (i in 1:(nf-1)){
	for (j in (i+1):nf){
		pairs[[length(pairs)+1]] <- c(i,j)
	}
}

for (p in pairs){
	pcx <- p[1]
        pcy <- p[2]
	
	plot <- ggplot(scores, aes(x=.data[[paste0("PC", pcx)]], y=.data[[paste0("PC", pcy)]], fill = SPECIES)) +
		geom_point(size=3, shape=21, color = "black") +
		geom_text_repel(aes(label=INDIVIDUALS), size=3, max.overlaps=Inf) +
		labs(x = paste0("PC", pcx, " (", percent_pca[paste0("PC",pcx)], "%)"), y = paste0("PC", pcy, " (", percent_pca[paste0("PC",pcy)], "%)"), title=paste0("PCA genetic structure for the ", method, " method"), fill = "Presumed species\nbased on field observations") +
		theme_classic() +
		theme(legend.position = "right", text=element_text(size=14), legend.title = element_text(size = 10), legend.text = element_text(size = 10))
	ggsave(paste0(outdir, "/", method, "_ACP_PC", pcx, "_PC", pcy, ".png"), plot)
	}

file.create(paste0(outdir, "/", method, "_ACP.done"))
