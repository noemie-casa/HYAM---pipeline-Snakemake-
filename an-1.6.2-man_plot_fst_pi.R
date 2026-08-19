#HYAM_stage
#Manhattan plot for all parameters

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

# Taille des chromosomes à partir du Fst
Fst$CHR <- factor(Fst$CHR, levels = sort(unique(Fst$CHR))) 

chr_sizes <- Fst %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POS))

chr_sizes <- chr_sizes %>%
  mutate(offset = c(0, cumsum(chr_len)[-n()]))

# Ajout des positions cumulées à chaque table
Fst <- Fst %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BPcum = POS + offset)

PiA <- PiA %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BPcum = POS + offset)

PiB <- PiB %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BPcum = POS + offset)

# Position des labels des chromosomes
axis_set <- Fst %>%
  group_by(CHR) %>%
  summarise(center = (max(BPcum) + min(BPcum)) / 2)

# couleurs 
chr_sizes$color_group <- rep(c("A", "B"), length.out = nrow(chr_sizes))
Fst <- Fst %>%
  left_join(chr_sizes[, c("CHR", "color_group")], by = "CHR")
PiA <- PiA %>%
  left_join(chr_sizes[, c("CHR", "color_group")], by = "CHR")
PiB <- PiB %>%
  left_join(chr_sizes[, c("CHR", "color_group")], by = "CHR")

# moyenne 
mean_Fst <- mean(Fst$VAL, na.rm=T)
mean_PiA <- mean(PiA$VAL, na.rm=T)
mean_PiB <- mean(PiB$VAL, na.rm=T)

# seuil 95%
q95_Fst <- quantile(Fst$VAL, probs = 0.95, na.rm=TRUE)
q95_PiA <- quantile(PiA$VAL, probs = 0.95, na.rm = TRUE)
q95_PiB <- quantile(PiB$VAL, probs = 0.95, na.rm = TRUE)

# plot
p1 <- ggplot(Fst, aes(BPcum, VAL, color=color_group)) +
	geom_point(size=1) +
	geom_hline(yintercept=mean_Fst, color="green", linewidth = 1, linetype = "solid") +
	geom_hline(yintercept=q95_Fst, color="red", linewidth = 1, linetype = "solid") +
	annotate("text", x=max(Fst$BPcum), y = max(Fst$VAL), label = paste0("Mean Fst = ", format(signif(mean_Fst, 3), scientific = T)), size = 6) +
	annotate("text", x=(max(Fst$BPcum)-0.05 * diff(range(Fst$BPcum))), y = (max(Fst$VAL)-0.05*diff(range(Fst$VAL))), label = paste0("Threshold quantile 95% = ", format(signif(q95_Fst, 3), scientific = T)), size = 6) +
	scale_color_manual(values=c("blue","orange")) +
	scale_x_continuous(labels=axis_set$CHR, breaks=axis_set$center) +
	theme_classic() +
	theme(legend.position="none") +
	labs(title=paste0("Fst between ", speciesA, " and ", speciesB), x=NULL, y="Fst")

p2 <- ggplot(PiA, aes(BPcum, VAL, color=color_group)) +
        geom_point(size=1) +
        geom_hline(yintercept=mean_PiA, color="green", linewidth = 1, linetype = "solid") +
        geom_hline(yintercept=q95_PiA, color="red", linewidth = 1, linetype = "solid") +
        annotate("text", x=max(PiA$BPcum), y = max(PiA$VAL), label = paste0("Mean Pi = ", format(signif(mean_PiA, 3), scientific = T)), size = 6) +
        annotate("text", x=(max(PiA$BPcum)-0.05 * diff(range(PiA$BPcum))), y = (max(PiA$VAL)-0.05*diff(range(PiA$VAL))), label = paste0("Threshold quantile 95% = ", format(signif(q95_PiA, 3), scientific = T)), size = 6) +
        scale_color_manual(values=c("blue","orange")) +
        scale_x_continuous(labels=axis_set$CHR, breaks=axis_set$center) +
        theme_classic() +
        theme(legend.position="none") +
        labs(title=paste0("Nucleotide diversity in ", speciesA), x=NULL, y="Pi")

p3 <- ggplot(PiB, aes(BPcum, VAL, color=color_group)) +
        geom_point(size=1) +
        geom_hline(yintercept=mean_PiB, color="green", linewidth = 1, linetype = "solid") +
        geom_hline(yintercept=q95_PiB, color="red", linewidth = 1, linetype = "solid") +
        annotate("text", x=max(PiB$BPcum), y = max(PiB$VAL), label = paste0("Mean Pi = ", format(signif(mean_PiB, 3), scientific = T)), size = 6) +
        annotate("text", x=(max(PiB$BPcum)-0.05 * diff(range(PiB$BPcum))), y = (max(PiB$VAL)-0.05*diff(range(PiB$VAL))), label = paste0("Threshold quantile 95% = ", format(signif(q95_PiB, 3), scientific = T)), size = 6) +
        scale_color_manual(values=c("blue","orange")) +
        scale_x_continuous(labels=axis_set$CHR, breaks=axis_set$center) +
        theme_classic() +
        theme(legend.position="none") +
        labs(title=paste0("Nucleotide diversity in ", speciesB), x=NULL, y="Pi")

png(output_file, width = 9000, height = 6000, res = 300)
(p1 / p2 / p3)
dev.off()

