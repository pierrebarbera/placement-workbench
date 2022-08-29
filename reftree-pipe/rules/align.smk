# =================================================================================================
#     Dummy for when input is already aligned
# =================================================================================================

rule no_alignment:
    group: "alignment"
    input:
        "{outdir}/result/{sample}/{autoref}/ref_candidates.fa"
    output:
        "{outdir}/result/{sample}/{autoref}/no_alignment/aligned.afa"
    log:
        "{outdir}/result/{sample}/{autoref}/no_alignment/log.txt"
    script:
        "../../common/symlink.py"
localrules: no_alignment

# =================================================================================================
#     Alignment with mafft
# =================================================================================================

rule align_mafft:
    group: "alignment"
    input:
        "{outdir}/result/{sample}/{autoref}/ref_candidates.fa"
    output:
        "{outdir}/result/{sample}/{autoref}/mafft/aligned.afa"
    params:
        auto        = True,
        maxiterate  = 1000
    threads:
        get_threads( "mafft" )
    log:
        "{outdir}/result/{sample}/{autoref}/mafft/alignment.log"
    conda:
        "../envs/mafft.yaml"
    script:
        "../scripts/mafft.py"

# =================================================================================================
#     Alignment with muscle
# =================================================================================================

rule align_muscle:
    group: "alignment"
    input:
        "{outdir}/result/{sample}/{autoref}/ref_candidates.fa"
    output:
        "{outdir}/result/{sample}/{autoref}/muscle/aligned.afa"
    params:
        extra = config["params"]["muscle"]["extra"]
    log:
        "{outdir}/result/{sample}/{autoref}/muscle/alignment.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        "muscle -in {input} -out {output} {params.extra} > {log} 2>&1"
