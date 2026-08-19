
#######################################
# Locus and SNP tables concatenation
#######################################
rule concatenation:
        input:
                file=expand(f"{Results}/DATA/{{sample}}_table.txt", sample=config["echantillons"].values())
        output:
                table=f"{Results}/DATA/{run_name}_table_completed.txt"
        params:
                outdir=f"{Results}/DATA"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                > {output.table}
                cat {input.file} >> {output.table}
                """


###############################
# Chromosome table creation 
###############################
rule CHROMcreation:
        input:
                fai=f"{ref_dir}/{fai}"
        output:
                table=f"{Results}/DATA/{species_ref}_chromosome_table.txt"
        params:
                outdir=f"{Results}/DATA"
        threads:
                1
        resources:
                mem_mb=64000
        run:
                chr = []
                with open(input.fai) as fai:
                        for line in fai:
                                fields = line.strip().split("\t")
                                chrom = fields[0]
                                size = int(fields[1])
                                chr.append((chrom, size))

                chr.sort(key=lambda x: x[1], reverse=True)

                with open(output.table, "w") as out:
                        for chrom, size in chr:
                                out.write(f"{chrom}\t1\t{size}\n")
                                


#################################
# Chromomap plot from R script 
#################################
rule ChromoMap:
        input:  
                table=f"{Results}/DATA/{run_name}_table_completed.txt",
                chromosome=f"{Results}/DATA/{species_ref}_chromosome_table.txt"
        output:
                map=f"{Results}/{run_name}_chrom_mapping.html"
        params:
                outdir=f"{Results}",
                script=f"{scriptdir}/ch-mapping.R",
                n_types=n_types,
                window_size=window_size
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # Module to load 
                module load r/4.5.2
                Rscript -e "library(chromoMap)"
                Rscript -e "library(htmlwidgets)"

                # Script execution
                Rscript {params.script} {input.table} {input.chromosome} {output.map} {params.n_types} {params.window_size}
                """
