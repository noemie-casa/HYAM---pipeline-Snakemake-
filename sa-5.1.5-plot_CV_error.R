#HYAM_stage
#Plot compared CV error

##################################################

#Packages to load
library(ggplot2)
library(tidyr)
library(dplyr)

#Files
args <- commandArgs(trailingOnly = TRUE)
error_file <- args[1]
output_file <- args[2]
method <- args[3]

cv <- read.table(error_file, header=F, sep=":", stringsAsFactors = F, col.names = c("K", "CV_error"))
cv$K <- as.numeric(gsub("K=", "", trimws(cv$K)))
cv$CV_error <- as.numeric(trimws(cv$CV_error))

#Plot
p <- ggplot(cv, aes(x=K, y=CV_error)) + 
	geom_line(linewidth = 1) +
	geom_point(size = 5) +
	scale_x_continuous(breaks = cv$K) +
	scale_y_continuous(limits = c(0, NA)) +
	labs(x = "Number of clusters (K)", y = "CV error", title = paste0("ADMIXTURE cross-validation error depending on number of K for the method ", method)) +
	theme_bw(base_size = 14) +
	theme(plot.title = element_text(size = 10))
ggsave(output_file,p) 
