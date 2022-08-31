# =================================================================================================
#     Cluster Query Sequences
# =================================================================================================

# special rule that skips clustering / does nothing
rule no_cluster:
    input:
        get_sample_fasta
    output:
        "{outdir}/no_clustering/samples/{sample}/queries.fa"
    log:
        "{outdir}/no_clustering/samples/{sample}/log.txt"
    script:
        "../../common/symlink.py"
localrules: no_cluster

rule cluster_swarm:
    group: "swarm"
    input:
        get_sample_fasta
    params:
        append_abundance = 1
    output:
        seeds           = "{outdir}/swarm/samples/{sample}/queries.fa",
        statistics_file = "{outdir}/swarm/samples/{sample}/otu_table.tsv"
    threads:
        get_threads( "swarm" )
    log:
        "{outdir}/swarm/samples/{sample}/log.txt"
    conda:
        "../envs/swarm.yaml"
    script:
        "../scripts/swarm.py"

rule cluster_dada2:
    input:
        get_sample_fastq
    output:
        fasta       = "{outdir}/dada2/samples/{sample}/queries.fa",
        otu_table   = "{outdir}/dada2/samples/{sample}/otu_table.tsv"
    log:
        "{outdir}/dada2/samples/{sample}/log.txt"
    conda:
        "../envs/dada2.yaml"
    script:
        "../scripts/dada2.R"
