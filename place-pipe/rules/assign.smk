# =================================================================================================
#     Perform Taxonomic Assignment
# =================================================================================================

rule assign:
    input:
        jplace      = "{outdir}/placed/{sample}.jplace",
        taxon_file  = config["data"]["taxonomy-file"]
    output:
        "{outdir}/assign/gappa/{sample}/profile.tsv"
    log:
        "{outdir}/logs/gappa/assign/{sample}.log"
    threads:
        get_highest_override( "gappa", "threads" )
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-assign.py"


# No need to execute this on the cluster computed nodes.
localrules: assign
