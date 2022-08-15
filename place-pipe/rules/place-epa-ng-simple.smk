# =================================================================================================
#     Setup Dependencies
# =================================================================================================

include: "place-epa-ng-common.smk"

# =================================================================================================
#     Placement with epa-ng
# =================================================================================================

rule epa_ng_place:
    input:
        tree = config["data"]["reference-tree"],
        msa  = config["data"]["reference-alignment"],
        sequences = "{outdir}/{clusterer}/aligned/{sample}.afa",

        # We here use a trick to allow the "model" to be either a RAxML model file, or a string.
        # Get the model param file from the config if it is a file, or, if that is empty,
        # request a file, which will be created by the raxml_ng_model_eval rule in epa-ng-common.
        # If it is neither a file nor empty, we assume it's an actual model string, which will
        # be used below in the params section then.
        model = epa_ng_place_model( "input" )
    output:
        jplace = protected("{outdir}/{clusterer}/placed/{sample}.jplace")
    params:
        # Get the model if it is a string (and not a file).
        model = epa_ng_place_model( "params" ),

        # Get any extra params to use with epa-ng
        extra = config["params"]["epa-ng"]["extra"]
    log:
        "{outdir}/{clusterer}/place/{sample}.log"
    conda:
        "../envs/epa-ng.yaml"
    threads:
        get_threads( "epa-ng" )
    script:
        "../scripts/epa-ng.py"
