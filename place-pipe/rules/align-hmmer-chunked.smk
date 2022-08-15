# =================================================================================================
#     Setup Dependencies
# =================================================================================================

include: "align-hmmer-common.smk"

# =================================================================================================
#     Alignment with hmmer
# =================================================================================================

rule hmmer_align:
    input:
        msa         = config["data"]["reference-alignment"],
        hmmprofile  = "{outdir}/hmmer/profile.hmm",
        sample      = "{outdir}/{clusterer}/chunkify/chunks/{sample}.fasta"
    output:
        sequences   = "{outdir}/{clusterer}/chunkify/aligned/{sample}.afa"
    params:
        extra       = config["params"]["hmmer"]["align-extra"],
        states      = hmmer_datatype_string
    log:
        "{outdir}/{clusterer}/aligned/{sample}.log"
    conda:
        "../envs/hmmer.yaml"
    threads:
        get_threads( "hmmer" )
    shell:
        "hmmalign"
        " --{params.states}"
        " --outformat afa"
        " -o {output.sequences}"
        " --mapali {input.msa}"
        " {params.extra} {input.hmmprofile} {input.sample}"
        " > {log} 2>&1"

