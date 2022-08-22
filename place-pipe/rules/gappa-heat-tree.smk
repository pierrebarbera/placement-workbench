# =================================================================================================
#     gappa heat-tree
# =================================================================================================

# Rule to create a heat tree visualization for all samples combined
rule gappa_heat_tree_all:
    input:
        expand( "{outdir}/{clusterer}/placed/{sample}.jplace",
                outdir=outdir,
                allow_missing=True,
                sample=sample_names
                )
    output:
        expand( "{outdir}/{clusterer}/heat-tree.{ext}",
                outdir=outdir,
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    params:
        extra=config["params"]["gappa"]["heat-tree"]["extra"],
    log:
        expand( "{outdir}/{clusterer}/heat_tree_all.{ext}.log",
                outdir=outdir,
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"
    group: "postplacement"

# Rule to create individual heat tree visualizations per sample
rule gappa_heat_tree:
    input:
        "{outdir}/placed/{sample}.jplace"
    output:
        expand( "{outdir}/{clusterer}/heat-trees/{{sample}}/heat-tree.{ext}",
                outdir=outdir,
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    params:
        extra=config["params"]["gappa"]["heat-tree"]["extra"],
    log:
        expand( "{outdir}/{clusterer}/heat_trees/{{sample}}.{ext}.log",
                outdir=outdir,
                allow_missing=True,
                ext=config["params"]["gappa"]["heat-tree"]["formats"]
                )
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"
    group: "postplacement"
