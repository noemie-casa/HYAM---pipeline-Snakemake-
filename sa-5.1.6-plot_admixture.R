#HYAM_stage
#Plot admixture

##################################################

#Packages to load
library(tidyverse)

#Files 
args <- commandArgs(trailingOnly = TRUE)
fam_file <- args[1]
Q_file <- args[2]
output_file <- args[3]
k <- as.numeric(args[4])

fam <- read.table(fam_file, sep="\t", header=F)
samples <- fam$V2

# load data
Qmat <- read.table(Q_file, header=F)
colnames(Qmat) <- paste0("Q", seq_len(k))
Qmat$sample <- samples
Qmat$k <- k

# order
Qcols <- paste0("Q", seq_len(k))
Qmat_ord <- Qmat %>% arrange(across(all_of(Qcols),desc))
sample_order <- Qmat_ord$sample

data <- Qmat %>% pivot_longer(cols = all_of(Qcols), names_to = "Q", values_to = "value") %>% mutate(sample = factor(sample, levels = sample_order))

# barplot
p <- ggplot(data, aes(x=sample, y=value, fill=factor(Q))) +
	geom_bar(stat="identity", position="stack") +
	xlab("Sample") +
	ylab("Ancestry proportion") +
	theme_bw() +
	theme(axis.text.x=element_text(angle=60, hjust=1))

if (k == 1){
	p <- p + scale_fill_brewer(palette="#66C2A5", name="Cluster")
} else {
	p <- p + scale_fill_brewer(palette="Set2", name="Cluster")
}

ggsave(filename=output_file, plot = p, height=8, width=24)
