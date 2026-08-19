#HYAM_stage
#Variant calling mapping

##################################################

#Packages to load
library(chromoMap)
library(htmlwidgets)

#Files 
args <- commandArgs(trailingOnly = TRUE)
SNP_file <- args[1]
chromosome_file <- args[2]
output_file <- args[3]
n_types <- as.numeric(args[4])
win_size <- as.numeric(args[5])

colors <- rainbow(n_types)

#Mapping
if( n_types == 1){
	SNP_mapping <- chromoMap(chromosome_file, SNP_file,
                         fixed.window=T,
                         window.size=win_size,
                         chr_color = "gainsboro",
                         legend = TRUE)

} else {
	SNP_mapping <- chromoMap(chromosome_file, SNP_file,
			 fixed.window=T,
			 window.size=win_size,
			 data_based_color_map = T,
			 data_type = "categorical", 
			 data_colors = list(c(colors)),
			 chr_color = "gainsboro",
                         legend = TRUE)
}
saveWidget(SNP_mapping, file = output_file)


