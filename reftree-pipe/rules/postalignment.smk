# =================================================================================================
#     Rules regarding post-alignment cleanup tasks
# =================================================================================================

# =================================================================================================
#     MSA trimming
# =================================================================================================

# trim ends that is always called before trimming occurs
rule clean_alignment:
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/aligned.afa"
    params:
        datatype    = config["settings"]["datatype"],
        n           = config["params"]["trim_ends_n"]
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/cleaned.afa"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/cleaner_log.txt"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/trim_ends.py"

# special rule that skips trimming / does nothing
rule no_trim:
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/cleaned.afa"
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/no_trim/trimmed.afa"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/no_trim/log.txt"
    script:
        "../../common/symlink.py"

rule trim_gblocks:
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/cleaned.afa"
    params:
        datatype    = ('p' if config["settings"]["datatype"] == 'aa' else 'd'),
        rel_input   = relative_input_path,
        extra       = config["params"]["gblocks"]["extra"]
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/gblocks/trimmed.afa"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/gblocks/log.txt"
    conda:
        "../envs/gblocks.yaml"
    shell:
        # somehow gblocks returns a non-zero exit value regardless of success or failure?!
        "$(Gblocks {input} -t={params.datatype} {params.extra} > {log} ; echo '' )"
        " && ln -s {params.rel_input}-gb {output}"

rule trim_trimal:
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/cleaned.afa"
    params:
        extra = config["params"]["trimal"]["extra"]
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/trimal/trimmed.afa"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/trimal/log.txt"
    conda:
        "../envs/trimal.yaml"
    shell:
        "trimal -in {input} -out {output} -fasta -automated1 {params.extra} 2> {log}"
localrules: no_trim, clean_alignment


# =================================================================================================
#     Remove duplicates
# =================================================================================================
