#HYAM_stage
#Biplot for Fst and Pi

##################################################

#Packages to load
library(ggplot2)
library(dplyr)
library(patchwork)

#Files 
args <- commandArgs(trailingOnly = TRUE)
Fst_file <- args[1]
PiA_file <- args[2]
PiB_file <- args[3]
output_file <- args[4]
speciesA <- args[5]
speciesB <- args[6]

Fst <- read.table(Fst_file, sep="\t", header=TRUE)
PiA <- read.table(PiA_file, sep="\t", header=TRUE)
PiB <- read.table(PiB_file, sep="\t", header=TRUE)

colnames(Fst) <- c("CHR", "POS", "VAL")
colnames(PiA) <- c("CHR", "POS", "VAL")
colnames(PiB) <- c("CHR", "POS", "VAL")

# moyenne
mean_Fst <- mean(Fst$VAL, na.rm=T)
mean_PiA <- mean(PiA$VAL, na.rm=T)
mean_PiB <- mean(PiB$VAL, na.rm=T)

# seuil 95%
q95_Fst <- quantile(Fst$VAL, probs = 0.95, na.rm=TRUE)
q5_PiA <- quantile(PiA$VAL, probs = 0.05, na.rm = TRUE)
q5_PiB <- quantile(PiB$VAL, probs = 0.05, na.rm = TRUE)

# dataframe
dataA <- inner_join(Fst, PiA, by = c("CHR", "POS"), suffix = c("_Fst", "_Pi"))
dataB <- inner_join(Fst, PiB, by = c("CHR", "POS"), suffix = c("_Fst", "_Pi"))

# labels légende
label_mean_Fst <- paste0("Mean Fst = ", format(signif(mean_Fst, 3), scientific = TRUE))
label_q95_Fst  <- paste0("Fst 95% = ", format(signif(q95_Fst, 3), scientific = TRUE))

label_mean_PiA <- paste0("Mean Pi = ", format(signif(mean_PiA, 3), scientific = TRUE))
label_q5_PiA   <- paste0("Pi 5% = ", format(signif(q5_PiA, 3), scientific = TRUE))

label_mean_PiB <- paste0("Mean Pi = ", format(signif(mean_PiB, 3), scientific = TRUE))
label_q5_PiB   <- paste0("Pi 5% = ", format(signif(q5_PiB, 3), scientific = TRUE))

# association labels et couleurs
colors_A <- c("green", "red", "green", "red")
names(colors_A) <- c(label_mean_PiA, label_q5_PiA, label_mean_Fst, label_q95_Fst)

colors_B <- c("green", "red", "green", "red")
names(colors_B) <- c(label_mean_PiB, label_q5_PiB, label_mean_Fst, label_q95_Fst)

# plot
pA <- ggplot(dataA, aes(x = VAL_Fst, y = VAL_Pi)) + 
	geom_point(size=1) +
        geom_hline(aes(yintercept=mean_PiA, color=label_mean_PiA), linewidth = 1, linetype = "solid") +
        geom_hline(aes(yintercept=q5_PiA, color=label_q5_PiA), linewidth = 1, linetype = "solid") +
	geom_vline(aes(xintercept=mean_Fst, color=label_mean_Fst), linewidth = 1, linetype = "dashed") +
        geom_vline(aes(xintercept=q95_Fst, color=label_q95_Fst), linewidth = 1, linetype = "dashed") +
        theme_classic() +
        theme(legend.position=c(0.90, 0.98), legend.justification = c(1, 1)) +
        labs(title=paste0("Nucleotide diversity of ", speciesA, " compared to the Fst between ", speciesA, " and ", speciesB), x="Fst", y="Pi") +
	scale_color_manual(values = colors_A, name = NULL)

pB <- ggplot(dataB, aes(x = VAL_Fst, y = VAL_Pi)) +
        geom_point(size=1) +
        geom_hline(aes(yintercept=mean_PiB, color=label_mean_PiB), linewidth = 1, linetype = "solid") +
        geom_hline(aes(yintercept=q5_PiB, color=label_q5_PiB), linewidth = 1, linetype = "solid") +
        geom_vline(aes(xintercept=mean_Fst, color=label_mean_Fst), linewidth = 1, linetype = "dashed") +
        geom_vline(aes(xintercept=q95_Fst, color=label_q95_Fst), linewidth = 1, linetype = "dashed") +
        theme_classic() +
        theme(legend.position=c(0.90, 0.98), legend.justification = c(1, 1)) +
        labs(title=paste0("Nucleotide diversity of ", speciesB, " compared to the Fst between ", speciesA, " and ", speciesB), x="Fst", y="Pi") +
        scale_color_manual(values = colors_B, name = NULL)


png(output_file, width = 9000, height = 6000, res = 300)
(pA / pB)
dev.off()

