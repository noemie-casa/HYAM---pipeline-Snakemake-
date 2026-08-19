#############################
# VCF for each population
#############################

###################################
# Bam list creation for freebayes
###################################
rule nuc_list_bam:
        input: 
                SAMPLES=lambda wildcards: expand(f"{bam_path}/{{sample}}_rdup.bam", sample = config["species"][wildcards.species])
        output:
                list=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.1-list_bam/{{species}}_bam_list.txt"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.1-list_bam"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                printf "%s\n" {input.SAMPLES} > {output.list}
                """


#####################################################################
# Variant calling - comparison between sample and reference genome
#####################################################################
rule nuc_comparaison:
        input:
                list=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.1-list_bam/{{species}}_bam_list.txt",
                REF=f"{ref_dir}/{genome}"
        output:
                vcf=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.2-variants/{{species}}_variants.vcf"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.2-variants",
                ploidy=ploidy,
                min_coverage=min_coverage
        threads:
                20
        resources:
                mem_mb=64000
        shell:
                """
                # module to load 
                module load freebayes/1.2.0
                module load samtools/1.9

                # indexation
                samtools faidx {input.REF}

                # comparison
                mkdir -p {params.outdir}
                freebayes \
                        --fasta-reference {input.REF} \
                        --ploidy {params.ploidy} \
                        --use-best-n-alleles 4 \
                        --min-coverage {params.min_coverage} \
                        --hwe-priors-off \
                        --bam-list {input.list} > {output.vcf}
                """


####################################################
# SNP calling - Selection of SNP among variants
####################################################
rule nuc_selection:
        input:
                variants=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.2-variants/{{species}}_variants.vcf"
        output:
                SNP=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.3-SNP/{{species}}_SNP.vcf"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.3-SNP",
                missing_data=missing_data
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load bcftools/1.16

                # selection
                mkdir -p {params.outdir}
                bcftools view -v snps -m 2 -M 2 -i 'F_MISSING < {params.missing_data} && MAF > 0' {input.variants} > {output.SNP}
                """


#######################################
# Nucleotide diversity calculation
#######################################
rule pi_windowed:
        input:
                vcf=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.3-SNP/{{species}}_SNP.vcf"
        output:
                pi_windowed=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{species}}_nuc_diversity.windowed.pi"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation",
                prefix=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{species}}_nuc_diversity",
                piwin=pi_win,
                pistep=pi_step
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load vcftools/0.1.16
                
                # nucleotide diversity
                mkdir -p {params.outdir}
                vcftools --vcf {input.vcf} --window-pi {params.piwin} --window-pi-step {params.pistep} --out {params.prefix}
                """


###########################################
# Nucletotide diversity table creation
###########################################
rule pi_table:
        input:
                pi=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{species}}_nuc_diversity.windowed.pi"
        output:
                table_center=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{species}}_nuc_diversity_table.txt"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi"
        threads:
                2
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk '
                BEGIN {{OFS="\t"}}
                {{mid=int(($2+$3)/2);print $1, mid, $5}} ' {input.pi} > {output.table_center}
                """


#####################################
# VCF files between two populations #
#####################################

###################################
# Bam list creation for freebayes
###################################
rule fst_list_bam:
        input: 
                SAMPLES=lambda wildcards: expand(f"{bam_path}/{{sample}}_rdup.bam", sample = (config["species"][wildcards.speciesA] + config["species"][wildcards.speciesB]))
        output:
                list=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.1-list_bam/{{speciesA}}_vs_{{speciesB}}_bam_list.txt"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.1-list_bam"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                printf "%s\n" {input.SAMPLES} > {output.list}
                """


#####################################################################
# Variant calling - comparison between sample and reference genome
#####################################################################
rule fst_comparaison:
        input:
                list=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.1-list_bam/{{speciesA}}_vs_{{speciesB}}_bam_list.txt",
                REF=f"{ref_dir}/{genome}"
        output:
                vcf=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.2-variants/{{speciesA}}_vs_{{speciesB}}_variants.vcf"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.2-variants",
                ploidy=ploidy,
                min_coverage=min_coverage
        threads:
                20
        resources:
                mem_mb=64000
        shell:
                """
                # module to load 
                module load freebayes/1.2.0
                module load samtools/1.9

                # indexation
                samtools faidx {input.REF}

                # comparison
                mkdir -p {params.outdir}
                freebayes \
                        --fasta-reference {input.REF} \
                        --ploidy {params.ploidy} \
                        --use-best-n-alleles 4 \
                        --min-coverage {params.min_coverage} \
                        --hwe-priors-off \
                        --bam-list {input.list} > {output.vcf}
                """


####################################################
# SNP calling - Selection of SNP among variants
####################################################
rule fst_selection:
        input:
                variants=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.2-variants/{{speciesA}}_vs_{{speciesB}}_variants.vcf"
        output:
                SNP=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.3-SNP/{{speciesA}}_vs_{{speciesB}}_SNP.vcf"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.3-SNP",
                missing_data=missing_data
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load bcftools/1.16

                # selection
                mkdir -p {params.outdir}
                bcftools view -v snps -m 2 -M 2 -i 'F_MISSING < {params.missing_data} && MAF > 0' {input.variants} > {output.SNP}
                """


##########################################
# Fixation index between two populations #
##########################################

##################################
# Indiv lists group by species
##################################
rule indiv_list:
        output:
                list=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list/{{species}}_indiv_list.txt"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list",
                individuals=lambda wildcards: " ".join(config["species"][wildcards.species])
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                for ind in {params.individuals} 
                do
                        printf "%s\n" "$ind" >> {output.list}
                done
                """

########################################
# Fixation index between two species
########################################
rule fixation_index_windowed:
        input:
                vcf=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.2-variants/{{speciesA}}_vs_{{speciesB}}_variants.vcf",
                list_speA=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list/{{speciesA}}_indiv_list.txt",
                list_speB=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list/{{speciesB}}_indiv_list.txt"
        output:
                Fst_window=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst.windowed.weir.fst"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst",
                prefix=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst",
                fstwin=fst_win,
                fststep=fst_step
        threads:
                2
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load vcftools/0.1.16

                # Execution
                mkdir -p {params.outdir}
                vcftools --vcf {input.vcf} --weir-fst-pop {input.list_speA} --weir-fst-pop {input.list_speB} --fst-window-size {params.fstwin} --fst-window-step {params.fststep} --out {params.prefix}
                """


#############################################
# Fst by window files cleaning (center)
#############################################
rule windowed_cleaning:
        input:
                fst=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst.windowed.weir.fst"
        output:
                fst_clean=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst/{{speciesA}}_vs_{{speciesB}}_clean_Fst_windowed.txt"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst"
        threads:
                2
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk '
                BEGIN {{
                        OFS="\t"
                }}
                {{
                        mid=int(($2+$3)/2);
                        print $1, mid, $5
                }}
                ' {input.fst} > {output.fst_clean}
                """


######################
# Analyses per genes #
######################
##################################
# Nucleotide diversity per SNP
##################################
rule pi:
        input:
                vcf=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.3-SNP/{{species}}_SNP.vcf"
        output:
                pi=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{species}}_nuc_diversity.sites.pi"
        params:
                outdir=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation",
                prefix=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{species}}_nuc_diversity"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load vcftools/0.1.16

                # nucleotide diversity
                mkdir -p {params.outdir}
                vcftools --vcf {input.vcf} --site-pi --out {params.prefix}
                """


##############################
# Fixation index per SNP
##############################
rule fixation_index:
        input:
                vcf=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.2-variants/{{speciesA}}_vs_{{speciesB}}_variants.vcf",
                list_speA=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list/{{speciesA}}_indiv_list.txt",
                list_speB=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.4-indiv_list/{{speciesB}}_indiv_list.txt"
        output:
                Fst=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst.weir.fst"
        params:
                outdir=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst",
                prefix=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst"
        threads:
                2
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load vcftools/0.1.16

                # Execution
                mkdir -p {params.outdir}
                vcftools --vcf {input.vcf} --weir-fst-pop {input.list_speA} --weir-fst-pop {input.list_speB} --out {params.prefix}
                """


#####################
# genes calling
#####################
rule genes:
        input: 
                anno=f"{ref_dir}/{annotation}"
        output:
                table=f"{Results}/2-analyses_by_genes/2.1-genes_tables/{species_ref}_genes_table.txt"
        params:
                outdir=f"{Results}/2-analyses_by_genes/2.1-genes_tables"
        threads:
                4
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk -F'\t' '
                BEGIN{{OFS="\t"}}

                # noms des chromosomes
                $3 == "region" {{
                        match($9, /Alias=([^;]+)/, arr)
                        alias[$1] = arr[1]
                        next
                }}

                # garder uniquement les gènes 
                $3 == "gene" {{
                        chr = ($1 in alias) ? alias[$1] : $1
                        match($9, /gene_id=([^;]+)/, gene)
                        print gene[1], chr, $4, $5
                }}
                ' {input.anno} > {output.table}
                """


##########################################
# association between SNPs and genes
##########################################
rule asso:
        input:
                snp_list=f"{SNPs_list}",
                gene_list=f"{Results}/2-analyses_by_genes/2.1-genes_tables/{species_ref}_genes_table.txt"
        output:
                table=f"{Results}/2-analyses_by_genes/2.1-genes_tables/snp_gene_table.txt"
        params:
                outdir=f"{Results}/2-analyses_by_genes/2.1-genes_tables"
        threads:
                4
        resources:
                mem_mb=64000
        run:
                import pandas as pd
           
                snp = pd.read_csv(input.snp_list, sep="\t", header=None, names=["CHROM", "POS"])
                snp["NAME"] = snp.index + 1
                gene = pd.read_csv(input.gene_list, sep="\t", header=None, names=["GENE", "CHROM", "BEG", "END"])

                merged = pd.merge(snp, gene, on="CHROM")
                selected = (merged["POS"] >= merged["BEG"]) & (merged["POS"] <= merged["END"])
                filter = merged[selected].copy()

                col = ["NAME","CHROM","POS","GENE"]

                exons = filter[col]

                exons.to_csv(output.table, sep="\t", index=False)


############################
# Fst and Pi tables
############################
rule merge_genes_tables:
        input:
                fst_file=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.5-Fst/{{speciesA}}_vs_{{speciesB}}_Fst.weir.fst",
                piA_file=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{speciesA}}_nuc_diversity.sites.pi",
                piB_file=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.4-pi_calculation/{{speciesB}}_nuc_diversity.sites.pi",
                snp_file=f"{Results}/2-analyses_by_genes/2.1-genes_tables/snp_gene_table.txt"
        output:
                tables=f"{Results}/2-analyses_by_genes/2.2-comparisons_tables/{{speciesA}}_vs_{{speciesB}}_fst_pi_genes_tables.csv"
        params:
                outdir=f"{Results}/2-analyses_by_genes/2.2-comparisons_tables"
        threads:
                8
        resources:
                mem_mb=64000
        run:
                import pandas as pd
                import numpy as np

                fst = pd.read_csv(input.fst_file, sep="\t", header=0, names=["CHROM", "POS", "FST"])
                piA = pd.read_csv(input.piA_file, sep="\t", header=0, names=["CHROM", "POS", f"PI_{wildcards.speciesA}"])
                piB = pd.read_csv(input.piB_file, sep="\t", header=0, names=["CHROM", "POS", f"PI_{wildcards.speciesB}"])

                snp = pd.read_csv(input.snp_file, sep="\t", header=0, names=["NAME", "CHROM", "POS", "GENE"])

                fst["FST"] = pd.to_numeric(fst["FST"], errors="coerce")
                piA[f"PI_{wildcards.speciesA}"] = pd.to_numeric(piA[f"PI_{wildcards.speciesA}"], errors="coerce")
                piB[f"PI_{wildcards.speciesB}"] = pd.to_numeric(piB[f"PI_{wildcards.speciesB}"], errors="coerce")

                metrics = (fst.merge(piA, on=["CHROM", "POS"], how="outer").merge(piB, on=["CHROM", "POS"], how="outer"))
                merged = pd.merge(snp, metrics, on=["CHROM", "POS"], how="inner")
                gene_stats = merged.groupby(["GENE", "CHROM"]).agg(
                    POS_START=('POS', 'min'),
                    POS_END=('POS', 'max'),
                    FST_MEAN=('FST', lambda x: x.mean(skipna=True)),
                    PIA_MEAN=(f"PI_{wildcards.speciesA}", lambda x: x.mean(skipna=True)),
                    PIB_MEAN=(f"PI_{wildcards.speciesB}", lambda x: x.mean(skipna=True))
                ).reset_index()
                gene_stats["POS_MID"] = ((gene_stats["POS_START"] + gene_stats["POS_END"]) / 2).astype(int)

                final_cols = ["GENE", "CHROM", "POS_MID", "FST_MEAN", "PIA_MEAN", "PIB_MEAN"]
                output_df = gene_stats[final_cols]
                output_df = output_df.rename(columns={"PIA_MEAN" : f"PI{wildcards.speciesA}_MEAN", "PIB_MEAN" : f"PI{wildcards.speciesB}_MEAN"}) 
                output_df.to_csv(output.tables, sep="\t", index=False)
                output_df = output_df.sort_values(["CHROM", "POS_MID"])
                output_df.to_csv(output.tables, sep="\t", index=False)

                
##########################################################
# Fst and Pi histogram distribution between two species 
###########################################################
rule fst_pi_hist:
        input:
                fst=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst/{{speciesA}}_vs_{{speciesB}}_clean_Fst_windowed.txt",
                piA=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesA}}_nuc_diversity_table.txt",
                piB=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesB}}_nuc_diversity_table.txt"
        output:
                distrib=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.1-histogram/{{speciesA}}_vs_{{speciesB}}_fst_pi_distrib_hist.png"
        params:
                outdir=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.1-histogram",
                script=f"{scriptdir}/an-1.6.1-hist_fst_pi.py"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load python/3.12

                # Script execution
                mkdir -p {params.outdir}
                python3 {params.script} {input.fst} {input.piA} {input.piB} {output.distrib} {wildcards.speciesA} {wildcards.speciesB}
                """


#####################################################
# Fst and Pi Manhattan plot between two species
#####################################################
rule fst_pi_man_plot:
        input:
                fst=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst/{{speciesA}}_vs_{{speciesB}}_clean_Fst_windowed.txt",
                piA=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesA}}_nuc_diversity_table.txt",
                piB=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesB}}_nuc_diversity_table.txt"
        output:
                plot=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.2-manhattan_plot/{{speciesA}}_vs_{{speciesB}}_fst_pi_man_plot.png"
        params:
                outdir=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.2-manhattan_plot",
                script=f"{scriptdir}/an-1.6.2-man_plot_fst_pi.R"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load r/4.5.2

                # Script execution
                mkdir -p {params.outdir}
                Rscript {params.script} {input.fst} {input.piA} {input.piB} {output.plot} {wildcards.speciesA} {wildcards.speciesB}
                """


############################################
# Pi and Fst tables between two species
############################################
rule merge_tables:
        input:
                fst_file=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst/{{speciesA}}_vs_{{speciesB}}_clean_Fst_windowed.txt",
                piA_file=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesA}}_nuc_diversity_table.txt",
                piB_file=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesB}}_nuc_diversity_table.txt"
        output:
                tables=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.3-table/{{speciesA}}_vs_{{speciesB}}_fst_pi_tables.csv"
        params:
                outdir=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.3-table"
        threads:
                8
        resources:
                mem_mb=64000
        run:
                import pandas as pd

                fst = pd.read_csv(input.fst_file, sep="\t", header=None, names=["CHROM", "POS", "FST"])
                piA = pd.read_csv(input.piA_file, sep="\t", header=None, names=["CHROM", "POS", f"PI_{wildcards.speciesA}"])
                piB = pd.read_csv(input.piB_file, sep="\t", header=None, names=["CHROM", "POS", f"PI_{wildcards.speciesB}"])

                piA[f"PI_{wildcards.speciesA}"] = pd.to_numeric(piA[f"PI_{wildcards.speciesA}"], errors="coerce")
                piB[f"PI_{wildcards.speciesB}"] = pd.to_numeric(piB[f"PI_{wildcards.speciesB}"], errors="coerce")

                merged = (fst.merge(piA, on=["CHROM", "POS"], how="outer").merge(piB, on=["CHROM", "POS"], how="outer"))
                merged["PI_MEAN"] = merged[[f"PI_{wildcards.speciesA}",f"PI_{wildcards.speciesB}"]].mean(axis=1, skipna=True)
                merged = merged.sort_values(["CHROM", "POS"])
                merged.to_csv(output.tables, sep="\t", index=False)


################################
# PiA vs Fst and PiB vs Fst 
################################
rule fst_pi_biplot:
        input:
                fst=f"{Results}/1-raw_analyse/1.5-fixation_index/1.5.6-clean_Fst/{{speciesA}}_vs_{{speciesB}}_clean_Fst_windowed.txt",
                piA=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesA}}_nuc_diversity_table.txt",
                piB=f"{Results}/1-raw_analyse/1.4-nucleotide_diversity/1.4.5-clean_pi/{{speciesB}}_nuc_diversity_table.txt"
        output:
                plot=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.4-biplot/{{speciesA}}_vs_{{speciesB}}_fst_pi_biplot.png"
        params:
                outdir=f"{Results}/1-raw_analyse/1.6-comparsion_2_species/1.6.4-biplot",
                script=f"{scriptdir}/an-1.6.4-biplot_fst_pi.R"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load r/4.5.2

                # Script execution
                mkdir -p {params.outdir}
                Rscript {params.script} {input.fst} {input.piA} {input.piB} {output.plot} {wildcards.speciesA} {wildcards.speciesB}
                """
