# =================================================================================================
#     Setup Dependencies
# =================================================================================================

include: "align-hmmer-common.smk"

# =================================================================================================
#     Alignment with hmmer
# =================================================================================================

rule hmmer_align:
    group: "alignment"
    input:
        msa      = config["data"]["reference-alignment"],
        hmmfile  = "{outdir}/hmmer/profile.hmm",
        seqfile  = "{outdir}/{clusterer}/chunkify/chunks/{sample}.fasta"
    output:
        "{outdir}/{clusterer}/chunkify/aligned/{sample}.afa"
    params:
        outformat   = "afa",
        states      = hmmer_datatype_string
    log:
        "{outdir}/{clusterer}/aligned/{sample}.log"
    conda:
        "../envs/hmmer.yaml"
    threads:
        get_threads( "hmmer" )
    script:
        "../scripts/hmmalign.py"
