#########################################
# Alignement on reference genome
#########################################
rule Alignment:
        input:
                unpack(get_captus),
                index=f"{ref_dir}/{genome}"
        output:
                bam=f"{Results}/1-alignment/{{captus}}.bam"
        params:
                outdir=f"{Results}/1-alignment"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load bwa-mem2/2.2.1 

                # alignement
                mkdir -p {params.outdir}
                bwa-mem2 mem -R "@RG\\tID:{wildcards.captus}\\tSM:{wildcards.captus}" {input.index} {input.fna} > {output.bam}
                """


#########################
# Bam file sorting
#########################
rule sorting:
        input:
                bam=f"{Results}/1-alignment/{{captus}}.bam"
        output:
                sorted_bam=f"{Results}/2-sorting/2.1-sorting/{{captus}}_sorted.bam"
        params:
                outdir=f"{Results}/2-sorting/2.1-sorting"
        threads:
                8
        resources:
                mem_mb=64000
        shell:
                """
                # module to load         
                module load samtools/1.9

                # sorting
                mkdir -p {params.outdir}
                samtools sort {input.bam} -o {output.sorted_bam}
                samtools index {output.sorted_bam}
                """


##############################
# Duplicats elimination
#############################
rule dupelimination:
        input:
                sorted_bam=f"{Results}/2-sorting/2.1-sorting/{{captus}}_sorted.bam"
        output:
                rdup=f"{Results}/2-sorting/2.2-rdup/{{captus}}_rdup.bam"
        params:
                outdir=f"{Results}/2-sorting/2.2-rdup"
        threads:
                8
        resources:
                mem_mb=64000
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


##########################
# Locus table creation
##########################
rule table:
        input:
                rdup=f"{Results}/2-sorting/2.2-rdup/{{captus}}_rdup.bam"
        output:
                table=f"{Results}/3-tables/{{captus}}_contigs_table.txt"
        params:
                outdir=f"{Results}/3-tables"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                # module to load
                module load samtools/1.9

                # table creation 
                mkdir -p {params.outdir}
                samtools view {input.rdup} | awk -v captus="{wildcards.captus}" '{{
			chr = $3;
                        start = $4;
                        cigar = $6;

                        sample = captus;
                        len = 0;
                        num = "";

                        for (i=1; i<=length(cigar); i++) {{
                                c = substr(cigar, i, 1);

                                if (c ~ /[0-9]/) {{
                                        num = num c;
                                }} else {{
                                        if (c ~ /[MDN=X]/) {{
                                                len += num;
                                        }}
                                        num = "";
                                }}
                        }}

                        end = start + len - 1;

                        locus = $1;

                        print locus "\t" chr "\t" start "\t" end "\t" sample;
                 }}' > {output.table}
                 """


###############################
# Locus tables concatenation
###############################
rule concatenation:
        input:
                sample=expand(f"{Results}/3-tables/{{captus}}_contigs_table.txt",captus=CAPTUS)
        output:
                table=f"{chromomap}/{ech}_contigs_table.txt"
        params:
                outdir=f"{chromomap}"
        threads:
                1
        resources:
                mem_mb=64000
        shell:
                """
                mkdir -p {params.outdir}
                cat {input.sample} >> {output.table}
                """

