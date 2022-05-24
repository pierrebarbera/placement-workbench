# =================================================================================================
#     Calculate Diversity Metrics
# =================================================================================================

rule diversity_guppy:
    input:
        "{outdir}/placed/{sample}.jplace"
    output:
        "{outdir}/diversity/guppy/{sample}.csv"
    params:
        include_pendant = ("--include-pendant" if config["params"]["guppy"]["fpd"]["include-pendant"] else ""),
        extra           = config["params"]["guppy"]["fpd"]["extra"]
    log:
        "{outdir}/logs/guppy/diversity/{sample}.log"
    conda:
        "../envs/pplacer.yaml"
    script:
        "guppy fpd -o {output} --csv {params.include_pendant} {input} > {log}"


# No need to execute this on the cluster computed nodes.
localrules: diversity_guppy
