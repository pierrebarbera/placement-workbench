rule rf_distances_between_samples:
    """
    Calculates RF-distances between all samples, taking the best tree from each
    """
    input:
        expand( "{outdir}/result/{sample}/best_result/raxml-ng/tree/best.newick",
                sample=sample_names,
                allow_missing=True
                )
    output:
        "{outdir}/result/intersample_rf_distances.txt"
    log:
        "{outdir}/result/rfdist.log"
    conda:
        "../envs/raxml-ng.yaml"
    script:
        "../scripts/raxml-ng-rfdist.py"

rule rf_distances_within_sample:
    """
    Calculates RF-distances between all best trees of a given sample (starting MSA/configuration)
    """
    input:
        expand( "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
                autoref=autoref_list,
                aligner=aligner_list,
                trimmer=trimmer_list,
                allow_missing=True
                )
    output:
        "{outdir}/result/{sample}/intrasample_rf_distances.txt"
    log:
        "{outdir}/result/{sample}/rfdist.log"
    conda:
        "../envs/raxml-ng.yaml"
    script:
        "../scripts/raxml-ng-rfdist.py"

rule rf_distances_within_treesearch:
    """
    Calculates RF-distances between all ml_trees of a given search / combination of tools
    """
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick"
    output:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/rf_distances.txt"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/rfdist.log"
    conda:
        "../envs/raxml-ng.yaml"
    script:
        "../scripts/raxml-ng-rfdist.py"
