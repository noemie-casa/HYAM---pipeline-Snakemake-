###########################
# CDS table creation
###########################
rule table :
        input: 
                anno=f"{ref_dir}/{annotation}"
        output:
                table=f"{chromomap}/{species_ref}_genes_table.txt"
        params:
                outdir=f"{chromomap}"
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
                        print gene, chr, $4, $5, $3
                }}
                ' {input.anno} > {output.table}
                """
