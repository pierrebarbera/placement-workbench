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
        sequences = "{outdir}/chunkify/aligned/{sample}.afa",

        # See the "simple" setup for an explanation of the model function that we use here.
        model = epa_ng_place_model( "input" )
    output:
        jplace = "{outdir}/chunkify/placed/{sample}.jplace"
    params:
        # Get the model if it is a string (and not a file).
        model = epa_ng_place_model( "params" ),

        # Get any extra params to use with epa-ng
        extra = config["params"]["epa-ng"]["extra"]
    log:
        "{outdir}/logs/place/{sample}.log"
    conda:
        "../envs/epa-ng.yaml"
    threads:
        get_highest_override( "epa-ng", "threads" )
    script:
        "../scripts/epa-ng.py"
