# =================================================================================================
#     gappa heat-tree
# =================================================================================================

# Rule to create a heat tree visualization for all samples combined
rule gappa_heat_tree_all:
    input:
        expand("{outdir}/placed/{sample}.jplace", outdir=outdir, sample=sample_names )
    output:
        expand("{outdir}/heat-tree.{ext}",
            outdir=outdir,
            ext=config["params"]["gappa"]["heat-tree"]["formats"]
            )
    params:
        extra=config["params"]["gappa"]["heat-tree"]["extra"],
    log:
        "{outdir}/logs/gappa/heat_tree_all.log"
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"

# Rule to create individual heat tree visualizations per sample
rule gappa_heat_tree:
    input:
        "{outdir}/placed/{sample}.jplace"
    output:
        expand("{outdir}/heat-trees/{{sample}}/heat-tree.{ext}",
            outdir=outdir,
            ext=config["params"]["gappa"]["heat-tree"]["formats"]
            )
    params:
        extra=config["params"]["gappa"]["heat-tree"]["extra"],
    log:
        "{outdir}/logs/gappa/heat_trees/{sample}.log"
    conda:
        "../envs/gappa.yaml"
    script:
        "../scripts/gappa-heat-tree.py"
