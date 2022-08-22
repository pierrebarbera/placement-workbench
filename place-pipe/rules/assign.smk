# =================================================================================================
#     Perform Taxonomic Assignment
# =================================================================================================

rule assign:
    input:
        jplace      = "{outdir}/{clusterer}/placed/{sample}.jplace",
        taxon_file  = config["data"]["taxonomy-file"]
    output:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/profile.tsv"
    log:
        "{outdir}/{clusterer}/taxonomic_assignment/gappa_assign/{sample}/log.txt"
    threads:
        get_threads( "gappa" )
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-assign.py"
    group: "postplacement"


# No need to execute this on the cluster computed nodes.
localrules: assign
