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
        seqfile  = "{outdir}/{clusterer}/samples/{sample}/queries.fa"
    output:
        "{outdir}/{clusterer}/aligned/{sample}.afa"
    params:
        outformat   = "afa",
        states      = hmmer_datatype_string
    log:
        "{outdir}/{clusterer}/aligned/{sample}.log"
    conda:
        "../envs/hmmer.yaml"
    script:
        "../scripts/hmmalign.py"