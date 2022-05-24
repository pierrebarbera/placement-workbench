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
        sample      = get_sample_fasta
    output:
        sequences   = "{outdir}/aligned/{sample}.afa"
    params:
        extra       = config["params"]["hmmer"]["align-extra"],
        states      = hmmer_datatype_string
    log:
        "{outdir}/logs/align/{sample}.log"
    conda:
        "../envs/hmmer.yaml"
    threads:
        get_highest_override( "hmmer", "threads" )
    shell:
        "hmmalign"
        " --{params.states}"
        " --outformat afa"
        " -o {output.sequences}"
        " --mapali {input.msa}"
        " {params.extra} {input.hmmprofile} {input.sample}"
        " > {log} 2>&1"
