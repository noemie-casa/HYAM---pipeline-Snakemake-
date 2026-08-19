#################
# Trimming
#################
rule trimming:
        input:
                R1=f"{DATA}/{{sample}}_R1.fastq.gz",
                R2=f"{DATA}/{{sample}}_R2.fastq.gz"
        output:
                R1=f"{Results}/0-trimming/{{sample}}_R1.clean.fastq",
                R2=f"{Results}/0-trimming/{{sample}}_R2.clean.fastq"
        params:
                outdir=f"{Results}/0-trimming",
                html=f"{Results}/0-trimming/{{sample}}_clean.html"
        threads:
                20
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load fastp/0.23.1
                
                # trimming
                mkdir -p {params.outdir}
                fastp -w 10 \
                        -i {input.R1} -I {input.R2} \
                        -o {output.R1} -O {output.R2} \
                        -h {params.html} \
                        -q 20 -u 20 -l 100 -y 30 -g 10 -x 10 -p 20
                """


#########################################
# Alignement on reference genome
#########################################
rule Alignment:
        input:
                R1=R1,
                R2=R2,
                index=f"{ref_dir}/{genome}"
        output:
                bam=temp(f"{Results}/1-alignment/{{sample}}.bam")
        params:
                outdir=f"{Results}/1-alignment"
        threads:
                20
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load bwa-mem2/2.2.1 

                # alignement
                mkdir -p {params.outdir}
                bwa-mem2 mem -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" {input.index} {input.R1} {input.R2} > {output.bam}
                """


######################
# Sorting by name
######################
rule name:
        input:
                bam=f"{Results}/1-alignment/{{sample}}.bam"
        output:
                bam_name=temp(f"{Results}/2-sorting/2.1-name/{{sample}}_name.bam")
        params: 
                outdir=f"{Results}/2-sorting/2.1-name"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                # module to load    
                module load samtools/1.9

                # conversion
                mkdir -p {params.outdir}
                samtools sort -n {input.bam} > {output.bam_name}
                """


############
# Fixmate
############
rule fixmate:
        input:
                bam_name=f"{Results}/2-sorting/2.1-name/{{sample}}_name.bam"
        output:
                bam_fix=temp(f"{Results}/2-sorting/2.2-fixmate/{{sample}}_fixmate.bam")
        params:
                outdir=f"{Results}/2-sorting/2.2-fixmate"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                # module to load    
                module load samtools/1.9

                # conversion
                mkdir -p {params.outdir}
                samtools fixmate -m {input.bam_name} {output.bam_fix} 
                """


############################
# Sorting by coordinate
############################
rule sorting:
        input:
                bam_fix=f"{Results}/2-sorting/2.2-fixmate/{{sample}}_fixmate.bam"
        output:
                sorted_bam=temp(f"{Results}/2-sorting/2.3-sorting/{{sample}}_sorted.bam")
        params:
                outdir=f"{Results}/2-sorting/2.3-sorting"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                # module to load         
                module load samtools/1.9

                # sorting
                mkdir -p {params.outdir}
                samtools sort {input.bam_fix} -o {output.sorted_bam}
                samtools index {output.sorted_bam}
                """


##############################
# Duplicats elimination
#############################
rule dupelimination:
        input:
                sorted_bam=f"{Results}/2-sorting/2.3-sorting/{{sample}}_sorted.bam"
        output:
                rdup=f"{Results}/2-sorting/2.4-rdup/{{sample}}_rdup.bam"
        params:
                outdir=f"{Results}/2-sorting/2.4-rdup"
        threads:
                1
        resources:
                mem_mb=128000
        shell:
                """
                # module to load
                module load samtools/1.9

                # elimination
                mkdir -p {params.outdir}
                samtools markdup -r {input.sorted_bam} {output.rdup}
                samtools index {output.rdup}
                rm -f *.bam
                """


###################################
# Bam list creation for freebayes
###################################
rule list:
        input: 
                SAMPLE=expand(f"{Results}/2-sorting/2.4-rdup/{{sample}}_rdup.bam",sample=SAMPLES)
        output:
                list=f"{Results}/2-sorting/2.5-list/all_indiv_bam_list.txt"
        params:
                outdir=f"{Results}/2-sorting/2.5-list"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                printf "%s\n" {input.SAMPLE} > {output.list}
                """


#####################################################################
# Variant calling - comparison between sample and reference genome
#####################################################################
rule comparaison:
        input:
                list=f"{Results}/2-sorting/2.5-list/all_indiv_bam_list.txt",
                REF=f"{ref_dir}/{genome}"
        output:
                vcf=f"{Results}/3-SNP_calling/3.1-variants/all_indiv_variants.vcf"
        params:
                outdir=f"{Results}/3-SNP_calling/3.1-variants",
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
rule SNP:
        input:
                variants=f"{Results}/3-SNP_calling/3.1-variants/all_indiv_variants.vcf"
        output:
                SNP=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        params:
                outdir=f"{Results}/3-SNP_calling/3.2-SNP",
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
                bcftools view -v snps -m 2 -M 2 -i 'F_MISSING < {params.missing_data}' {input.variants} > {output.SNP}
                """


#######################################
# Individuals missing data
#######################################
rule ind_missing_data:
        input:
                SNP=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                missing=f"{Results}/3-SNP_calling/3.3-missing_data/missing_data_per_individual.imiss"
        params:
                outdir=f"{Results}/3-SNP_calling/3.3-missing_data"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load vcftools/0.1.16

                mkdir -p {params.outdir}
                vcftools --vcf {input.SNP} --missing-indv --out tmp_missing
                mv tmp_missing.imiss {output.missing}
                """


###########################
# SNP missing data
###########################
rule SNP_missing_data:
        input:
                SNP=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                missing=f"{Results}/3-SNP_calling/3.3-missing_data/missing_data_per_SNP.lmiss"
        params:
                outdir=f"{Results}/3-SNP_calling/3.3-missing_data"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load vcftools/0.1.16

                mkdir -p {params.outdir}
                vcftools --vcf {input.SNP} --missing-site --out tmp_missing
                mv tmp_missing.lmiss {output.missing}
                """


############################
# Plot SNP missing data
############################
rule hist_missingSNPdata:
        input:
                missing=f"{Results}/3-SNP_calling/3.3-missing_data/missing_data_per_SNP.lmiss"
        output:
                plot=f"{Results}/3-SNP_calling/3.3-missing_data/plot_missing_data_per_SNP.png"
        params:
                outdir=f"{Results}/3-SNP_calling/3.3-missing_data",
                script=f"{scriptdir}/sa-3.3-hist_SNP_missingdata.py"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load python/3.12

                # script execution
                mkdir -p {params.outdir}
                python3 {params.script} {input.missing} {output.plot}
                """


###########################
# SNP table creation 
###########################
rule SNPcreation:
        input:
                vcf=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                table=f"{chromomap}/all_indiv_SNP_table.txt"
        params:
                outdir=f"{chromomap}"
        threads:
                4
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk '
                BEGIN {{ OFS="\t"; c=0 }}
                /^#/ {{ next }}
                {{c++; print "SNP_"c, $1, $2, $2, "SNP"}}
                ' {input.vcf} > {output.table}
                """


################
# SNPs list
################
rule raw_SNP_list:
        input:
                vcf=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                list=f"{Results}/4-thinning/4.1-raw_analyse/4.1.1-SNP_list/no_filtered_SNP_all_species.txt"
        params:
                outdir=f"{Results}/4-thinning/4.1-raw_analyse/4.1.1-SNP_list"
        threads:
                4
        resources:
                mem_mb=64000
        run:
                records = []
                with open(input.vcf, "r") as vcf:
                        for line in vcf:
                                if line.startswith("#"):
                                        continue
                                fields = line.strip().split("\t")
                                chrom = fields[0]
                                pos = int(fields[1])
                                records.append((chrom, pos))
                
                with open(output.list, "w") as out:
                         for chrom, pos in records:
                                out.write(f"{chrom}\t{pos}\n")



##################################
# Raw analyses and visualisation #
##################################

######################################
# Distance matrix from python script 
######################################
rule raw_dist_mat:
        input:
                list=f"{Results}/4-thinning/4.1-raw_analyse/4.1.1-SNP_list/no_filtered_SNP_all_species.txt"
        output:
                distance=f"{Results}/4-thinning/4.1-raw_analyse/4.1.2-distance/all_SNP_dist_mat.txt"
        params:
                outdir=f"{Results}/4-thinning/4.1-raw_analyse/4.1.2-distance",
                script=f"{scriptdir}/an-1.2-comp_dist_mat.py"
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
                python3 {params.script} {input.list} {output.distance}
                """


######################################
# Distance distribution histogram 
######################################
rule raw_hist_dist:
        input:
                distance=f"{Results}/4-thinning/4.1-raw_analyse/4.1.2-distance/all_SNP_dist_mat.txt"
        output:
                distrib=f"{Results}/4-thinning/4.1-raw_analyse/4.1.2-distance/all_SNP_dist_hist.png"
        params:
                outdir=f"{Results}/4-thinning/4.1-raw_analyse/4.1.2-distance",
                script=f"{scriptdir}/an-1.2-hist_dist.py",
                method="no selective"
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
                python3 {params.script} {input.distance} {output.distrib} "{params.method}"
                """


######################################
# Distance distribution histogram 
######################################
rule raw_nbr_SNPs:
        input:
                list=f"{Results}/4-thinning/4.1-raw_analyse/4.1.1-SNP_list/no_filtered_SNP_all_species.txt"
        output:
                tab=f"{Results}/4-thinning/4.1-raw_analyse/4.1.3-nbr_SNP/nbr_SNPs_chromosome_all_SNP.tsv"
        threads:
                1
        resources:
                mem_mb=64000
        run:
                # Module to load 
                import pandas as pd
                import os
                all_counts = []

                # Compte des SNPs pour chaque chromosome
                file = input.list
                df = pd.read_csv(file, sep="\t", header=None, names=["CHROM", "POS"])
                counts = (df.groupby("CHROM").size().reset_index(name="NBR_SNP"))
                all_counts.append(counts)

                # Génération du fichier de sortie
                all_counts = pd.concat(all_counts)
                all_counts.to_csv(output.tab, sep = "\t", index=False)


############
# thinning #
############
rule no_chrom_position:
        input:
                vcf=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf",
                fai=f"{ref_dir}/{fai}"
        output:
                cum_pos=f"{Results}/4-thinning/cumulative_positions.txt"
        params:
                outdir=f"{Results}/4-thinning"
        threads:
                8
        resources:
                mem_mb=64000
        run:
                # chromosome size
                chrom_sizes = {}
                offsets = {}
                with open(input.fai) as fai:
                        for line in fai:
                                fields = line.strip().split("\t")
                                chrom = fields[0]
                                size = int(fields[1])
                                offset = int(fields[2])
                                chrom_sizes[chrom] = size
                                offsets[chrom] = offset

                # lecture du vcf
                records = []
                with open(input.vcf, "r") as vcf:
                        for line in vcf:
                                if line.startswith("#"):
                                        continue
                                fields = line.strip().split("\t")
                                chrom = fields[0]
                                pos = int(fields[1])
                                pos_cum = offsets[chrom] + pos
                                records.append((chrom, pos, pos_cum))


                # tri
                records.sort(key=lambda x: x[2])

                # écriture du fichier de sortie
                with open(output.cum_pos, "w") as out:
                        out.write("CHROM\tPOS\tPOS_CUM\n")
                        for chrom, pos, pos_cum in records:
                                out.write(f"{chrom}\t{pos}\t{pos_cum}\n")


#######################
# random thinning
#######################
rule random_SNP_selection:
        input: 
                cum_pos=f"{Results}/4-thinning/cumulative_positions.txt"
        output:
                filter=f"{Results}/4-thinning/random/filtered_SNP_random.txt"
        params:
                outdir=f"{Results}/4-thinning/random",
                nsnp=n_SNPs,
                dist_min=gap_SNP
        threads:
                8
        resources:
                mem_mb=64000
        run:
                # module to import
                import random as rd

                # lecture des SNPs
                records = []
                with open(input.cum_pos) as f:
                        next(f)
                        for line in f :
                                chrom, pos, pos_cum = line.strip().split("\t")
                                records.append((chrom, int(pos), int(pos_cum)))

                # distance moyenne
                distances = []
                for i in range (1, len(records)):
                        previous_pos = records[i-1][2]
                        current_pos = records[i][2]
                        distances.append(current_pos - previous_pos)
                dist_moy = sum(distances)/len(distances)

                # taille des fenêtres
                d = max(1, len(records)//params.nsnp)
                window_size = max(1, int(d * dist_moy))

                # Sélection des SNPs
                windows = {}
                for rec in records:
                        chrom, pos, pos_cum = rec
                        win = int(pos_cum // window_size)
                        if win not in windows:
                                windows[win] = []
                        windows[win].append(rec)
                
                selected = []
                previous_pos_cum = None
                for win in sorted(windows):
                        candidates = windows[win].copy()
                        rd.shuffle(candidates)
                        for snp in candidates:
                                current_pos_cum = snp[2]
                                if (previous_pos_cum is None or abs(current_pos_cum - previous_pos_cum) >= params.dist_min):
                                        selected.append(snp)
                                        previous_pos_cum = current_pos_cum

                                        break

                # Tri
                selected.sort(key=lambda x : x[2])

                # Ecriture du fichier de sortie
                with open(output.filter, "w") as out:
                        for chrom, pos, pos_cum in selected:
                                out.write(f"{chrom}\t{pos}\n")


#############################
# Proportional thinning
#############################
rule proportional_SNP_selection:
        input: 
                cum_pos=f"{Results}/4-thinning/cumulative_positions.txt"
        output:
                filter=f"{Results}/4-thinning/proportional/filtered_SNP_proportional.txt"
        params:
                outdir=f"{Results}/4-thinning/proportional",
                nsnp=n_SNPs,
                dist_min=gap_SNP
        threads:
                8
        resources:
                mem_mb=64000
        run:
                # module to import
                import random as rd
                from collections import defaultdict

                # lecture des SNPs
                records = []
                with open(input.cum_pos) as f:
                        next(f)
                        for line in f :
                                chrom, pos, pos_cum = line.strip().split("\t")
                                records.append((chrom, int(pos), int(pos_cum)))

                # distance moyenne
                distances = []
                for i in range (1, len(records)):
                        previous_pos = records[i-1][2]
                        current_pos = records[i][2]
                        distances.append(current_pos - previous_pos)
                dist_moy = sum(distances)/len(distances)

                # taille des fenêtres
                d = max(1, len(records)//params.nsnp)
                window_size = max(1, int(d * dist_moy))

                # Sélection des SNPs
                windows = {}
                for rec in records:
                        chrom, pos, pos_cum = rec
                        win = int(pos_cum // window_size)
                        if win not in windows:
                                windows[win] = []
                        windows[win].append(rec)
                
                selected = []
                for win in sorted(windows):
                        candidates = windows[win].copy()
                        density = len(candidates) / len(records)
                        n_select = round(density * params.nsnp)
                        if n_select == 0:
                                continue
                        rd.shuffle(candidates)
                        local_count = 0

                        for snp in candidates:
                                current_pos = snp[2]
                                too_close = False
                                for kept in selected:
                                        if abs(current_pos - kept[2]) < params.dist_min:
                                                too_close = True
                                                break
                                if not too_close:
                                        selected.append(snp)
                                        local_count += 1
                                if local_count >= n_select:
                                        break
                        
                # Tri
                selected.sort(key=lambda x : x[2])

                # Ecriture du fichier de sortie

                with open(output.filter, "w") as out:
                        for chrom, pos, pos_cum in selected:
                                out.write(f"{chrom}\t{pos}\n")


##################
# ACP thinning
##################
rule ACP_SNP_selection:
        input:
                vcf=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                filter=f"{Results}/4-thinning/ACP/filtered_SNP_ACP.txt"
        params:
                outdir=f"{Results}/4-thinning/ACP",
                nsnp=n_SNPs,
                dist_min=gap_SNP,
                script=f"{scriptdir}/sa-4-thinning_ACP.R",
                n_axes=n_axes
        threads:
                8
        resources:
                mem_mb=64000
        shell: 
                """
                # Module to load
                module load r/4.5.2
                Rscript -e "library(adegenet)"
                Rscript -e "library(vcfR)"

                # Script execution
                mkdir -p {params.outdir}
                Rscript {params.script} {input.vcf} {output.filter} {params.nsnp} {params.dist_min} {params.n_axes}
                """


rule SNP_table:
        input:
                list=f"{Results}/4-thinning/{{method}}/filtered_SNP_{{method}}.txt"
        output:
                table=f"{chromomap}/{{method}}_filtered_SNP_table.txt"
        params:
                outdir=f"{chromomap}"
        threads:
                4
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk '
                BEGIN {{ OFS="\t"; c=0 }}

                /^#/ {{ next }}
                {{c++; print "SNP_"c, $1, $2, $2, "SNP"}}
                ' {input.list} > {output.table}
                """

rule bgzip_vcf:
        input:
                vcf=f"{Results}/3-SNP_calling/3.2-SNP/all_indiv_SNP.vcf"
        output:
                vcf=f"{Results}/4-thinning/{{method}}/index_trio.vcf.gz",
                tbi=f"{Results}/4-thinning/{{method}}/index_trio.vcf.gz.tbi"
        params:
                outdir=f"{Results}/4-thinning/{{method}}"
        threads:
                4
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load
                module load bcftools/1.16

                mkdir -p {params.outdir}
                bgzip -c {input.vcf} > {output.vcf}
                bcftools index -t {output.vcf}
                """

rule filtered_vcf:
        input:
                vcf=f"{Results}/4-thinning/{{method}}/index_trio.vcf.gz",
                tbi=f"{Results}/4-thinning/{{method}}/index_trio.vcf.gz.tbi",
                snp=f"{Results}/4-thinning/{{method}}/filtered_SNP_{{method}}.txt"
        output:
                vcf=f"{Results}/4-thinning/{{method}}/{{method}}_selected_snps.vcf"
        params:
                outdir=f"{Results}/4-thinning/{{method}}"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load 
                module load bcftools/1.16
               
                # file creation 
                bcftools view -R {input.snp} {input.vcf} -o {output.vcf}
                """
                
                
#############
# ADMIXTURE #
#############

######################################
# Input compilation for ADMIXTURE
######################################
rule chrom_list:
        input:
                fai=f"{ref_dir}/{fai}"
        output:
                list=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list/chrom_list.txt"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list"
        threads:
                1
        resources:
                mem_mb=64000
        run:
                chrom = []

                with open(input.fai) as fai:
                        for line in fai:
                                fields = line.strip().split("\t")
                                chrom.append(fields[0])

                with open(output.list, "w") as out:
                        for c in chrom :
                                 out.write(f"{c}\n")

rule chrom_corres:
        input:
                ref=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list/chrom_list.txt"
        output:
                corres=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list/chrom_correspondance.txt"
        params: 
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                awk '{{print $0 "\t1"}}' {input.ref} > {output.corres}
                """

rule chrom_anno:
        input:
                vcf=f"{Results}/4-thinning/{{method}}/{{method}}_selected_snps.vcf",
                corres=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.1-chrom_list/chrom_correspondance.txt"
        output:
                anno=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.2-vcf_anno_chrom/renamed_chr_{{method}}_SNP.vcf"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.2-vcf_anno_chrom"
        threads:
                4
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load
                module load bcftools/1.16
  
                # Annotation
                mkdir -p {params.outdir}
                bcftools annotate --rename-chrs {input.corres} {input.vcf} > {output.anno}
                """

rule vcf2plink:
        input:
                vcf=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.2-vcf_anno_chrom/renamed_chr_{{method}}_SNP.vcf"
        output:
                bed=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP.bed",
                bim=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP.bim",
                fam=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP.fam"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink",
                prefixe=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load
                module load plink2/2.00a5.12
 
                # Convertion
                mkdir -p {params.outdir}
                plink2 --allow-extra-chr --chr-set 1 --make-bed --vcf {input.vcf} --out {params.prefixe}

                """


rule admixture:
        input:
                bed=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP.bed"
        output:
                log=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results/{{method}}_log{{k}}.out",
                P=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results/renamed_chr_{{method}}_SNP.{{k}}.P",
                Q=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results/renamed_chr_{{method}}_SNP.{{k}}.Q"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load
                module load admixture/1.3.0

                # log cration
                mkdir -p {params.outdir}
                admixture --cv {input.bed} {wildcards.k} > {output.log}
                mv renamed_chr_{wildcards.method}_SNP.{wildcards.k}.P {output.P}
                mv renamed_chr_{wildcards.method}_SNP.{wildcards.k}.Q {output.Q}
                """

rule k_best_choice:
        input:
                expand(f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results/{{method}}_log{{k}}.out", k=range(kmin, (kmax+1)), method=METHOD)
        output:
                error=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.5-CV_error/{{method}}_renamed_chr_SNP.cv.error"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.5-CV_error",
                inputdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                grep "CV error" {params.inputdir}/{wildcards.method}_log*.out | awk -F'[=)]' '{{print "K="$2, $NF}}' > {output.error}
                """

rule plot_cv_error:
        input:
                error=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.5-CV_error/{{method}}_renamed_chr_SNP.cv.error"
        output:
                plot=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.5-CV_error/plot_{{method}}_renamed_chr_SNP_cv_error.png"
        params: 
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.5-CV_error",
                script=f"{scriptdir}/sa-5.1.5-plot_CV_error.R"
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
                Rscript {params.script} {input.error} {output.plot} {wildcards.method}
                """


#######################
# Plot ADMIXTURE
#######################
rule admix_plot:
        input:
                Q=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.4-results/renamed_chr_{{method}}_SNP.{{k}}.Q",
                fam=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.3-plink/renamed_chr_{{method}}_SNP.fam"
        output:
                plot=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.6-plot/{{method}}_admixture_plot_k{{k}}.png"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.1-ADMIXTURE/5.1.6-plot",
                script=f"{scriptdir}/sa-5.1.6-plot_admixture.R"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load r/4.5.2
                Rscript -e "library(tidyverse)"

                # Script execution
                mkdir -p {params.outdir}
                Rscript {params.script} {input.fam} {input.Q} {output.plot} {wildcards.k}
                """


#########################################################
# Presumed species information for each individual
#########################################################
rule presumed_species:
        output:
                list=f"{Results}/5-after_thinning_analyse/5.2-ACP/5.2.1-species_list/presumed_species_information.txt"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.2-ACP/5.2.1-species_list"
        threads:
                1
        resources:
                mem_mb=64000
        run:
                import os
                os.makedirs(params.outdir, exist_ok=True)

                with open(output.list, "w") as out:
                        out.write("INDIVIDUALS\tSPECIES\n")

                        for species, individuals in config["species"].items():
                                for ind in individuals:
                                        out.write(f"{ind}\t{species}\n")


##################
# Plot ACP
##################
rule ACP_plot:
        input:
                vcf=f"{Results}/4-thinning/{{method}}/{{method}}_selected_snps.vcf",
                species=f"{Results}/5-after_thinning_analyse/5.2-ACP/5.2.1-species_list/presumed_species_information.txt"
        output:
                touch(f"{Results}/5-after_thinning_analyse/5.2-ACP/5.2.2-ACP_plot/{{method}}_ACP.done")
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.2-ACP/5.2.2-ACP_plot",
                script=f"{scriptdir}/sa-5.2.2-ACP.R",
                n_axes=n_axes,
                method=lambda wildcards: wildcards.method
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load r/4.5.2
                Rscript -e "library(ggplot2)"
                Rscript -e "library(adegenet)"
                Rscript -e "library(ggrepel)"
                Rscript -e "library(vcfR)"

                # Script execution
                mkdir -p {params.outdir}
                Rscript {params.script} {input.vcf} {params.outdir} {params.n_axes} {params.method} {input.species}
                """


#########################################
# Comparaison between different methods #
#########################################

######################################
# Distance matrix from pyhton script 
######################################
rule comp_dist_mat:
        input:
                filtered=f"{Results}/4-thinning/{{method}}/filtered_SNP_{{method}}.txt"
        output:
                distance=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.1-file_generation/{{method}}_dist_mat.txt"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.1-file_generation",
                script=f"{scriptdir}/an-1.2-comp_dist_mat.py"
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
                python3 {params.script} {input.filtered} {output.distance}
                """


######################################
# Distance distribution histogram 
######################################
rule comp_hist_dist:
        input:
                distance=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.1-file_generation/{{method}}_dist_mat.txt"
        output:
                distrib=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.2-histogramme_distance_SNPs/{{method}}_dist_hist.png"
        params:
                outdir=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.2-histogramme_distance_SNPs",
                script=f"{scriptdir}/an-1.2-hist_dist.py"
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
                python3 {params.script} {input.distance} {output.distrib} {wildcards.method}
                """

#################################################
# Number of SNPs per chromosome for each method
#################################################
rule comp_nbr_SNPs:
        input:
                filtered=expand(f"{Results}/4-thinning/{{method}}/filtered_SNP_{{method}}.txt", method=METHOD)
        output:
                tab=f"{Results}/5-after_thinning_analyse/5.3-comparaison/5.3.1-file_generation/nbr_SNPs_chromosome_methode.tsv"
        threads:
                1
        resources:
                mem_mb=64000
        run:
                # Module to load 
                import pandas as pd
                import os

                all_counts = []

                # Compte des SNPs pour chaque chromosome pour chaque méthode
                for file in input.filtered:
                        method = os.path.basename(file)
                        method = method.replace("filtered_SNP_", "")
                        method = method.replace(".txt", "")
                        df = pd.read_csv(file, sep="\t", header=None, names=["CHROM", "POS"])
                        counts = (df.groupby("CHROM").size().reset_index(name="NBR_SNP"))
                        counts["METHOD"] = method
                        all_counts.append(counts)

                # Génération du fichier de sortie
                all_counts = pd.concat(all_counts)
                result = all_counts.pivot(index="CHROM", columns="METHOD", values="NBR_SNP").fillna(0).astype(int)
                result.to_csv(output.tab, sep = "\t")
