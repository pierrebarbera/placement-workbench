# =================================================================================================
#     Calculate Diversity Metrics
# =================================================================================================

rule diversity_guppy:
    group: "postplacement"
    input:
        "{outdir}/{clusterer}/placed/{sample}.jplace"
    output:
        "{outdir}/{clusterer}/diversity/guppy_fpd/{sample}.csv"
    params:
        include_pendant = ("--include-pendant" if config["params"]["guppy"]["fpd"]["include-pendant"] else ""),
        extra           = config["params"]["guppy"]["fpd"]["extra"]
    log:
        "{outdir}/{clusterer}/diversity/guppy_fpd/{sample}.log"
    conda:
        "../envs/pplacer.yaml"
    shell:
        "guppy fpd -o {output} --csv {params.include_pendant} {input} > {log}"

# # No need to execute this on the cluster computed nodes.
# localrules: diversity_guppy
