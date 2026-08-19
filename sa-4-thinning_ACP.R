#HYAM_stage
#ACP_thining

##################################################

#Packages to load
library(adegenet)
library(vcfR)

# Files 
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- args[2]
n_SNP <- as.numeric(args[3])
dist_min <- as.numeric(args[4])
nf <- as.numeric(args[5])

# Import vcf
vcf <- read.vcfR(vcf_file)
genlight <- vcfR2genlight(vcf)

# PCA et poids de chaque axe
pca_res <- glPca(genlight, nf = nf, loadings = TRUE)
percent_pca <- pca_res$eig / sum(pca_res$eig)

# SNP loadings
loadings <- pca_res$loadings[, 1:nf]

# Scores pondérés
score <- rowSums(abs(loadings) * matrix(percent_pca[1:nf], nrow=nrow(loadings), ncol=nf, byrow=T))
score_df <- data.frame(SNP = rownames(loadings), SCORE=score)

# Positions des SNPs
snp_pos <- data.frame(SNP=rownames(loadings), CHROM=vcf@fix[,"CHROM"], POS=as.numeric(vcf@fix[,"POS"]))

# Tri
score_df <- merge(score_df, snp_pos, by = "SNP")
score_df <- score_df[order(score_df$SCORE, decreasing = TRUE), ]

# Sélection
selected_snps <- data.frame()
for(i in seq_len(nrow(score_df))){
	snp <- score_df[i, ]
	if (nrow(selected_snps) == 0) {
		selected_snps <- rbind(selected_snps, snp)
	} else {
		same_chr <- selected_snps$CHROM == snp$CHROM
		dist_ok <- (sum(same_chr) == 0) || all(abs(selected_snps$POS[same_chr] - snp$POS) >= dist_min)
		if (dist_ok) {
			selected_snps <- rbind(selected_snps, snp)
		}
	}
	if (nrow(selected_snps) >= n_SNP){
		break
	}
}

final_out <- selected_snps[, c("CHROM", "POS")]
final_out <- final_out[order(final_out$CHROM, final_out$POS), ]
write.table(final_out, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
