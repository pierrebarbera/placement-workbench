# =================================================================================================
#     Auxilliary Rules and Functions
# =================================================================================================

# Rule to infer the best ML model parameters given the MSA and tree that we are going to use
# for placement.
rule raxml_ng_model_eval:
    input:
        tree = config["data"]["reference-tree"],
        msa = config["data"]["reference-alignment"]
    output:
        modelfile = "{outdir}/model/model_eval.raxml.bestModel",
        model_dir = directory( "{outdir}/model" )
    params:
        model = config["params"]["epa-ng"]["model"]
    log:
        "{outdir}/model/model_eval.log"
    threads:
        get_threads( "raxml-ng" )
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --evaluate"
        " --msa {input.msa}"
        " --tree {input.tree}"
        " --model {params.model}"
        " --prefix {output.model_dir}/model_eval"
        " --threads {threads}"
        " > {log} 2>&1"
    group: "placement"

# epa-ng expects a model of evolution to be given either as a string (such as the model
# string used in raxml-ng), or as a file (such as the best model output file of raxml-ng).
# For snakemake, this is tricky, because if we are given a file, we want this to exist.
# If its a string, not. If it's emtpy however, we create the file ourself first (via the above rule),
# and so also want to point to a file.
# So let's use a function to make this distinction:
# The function is called twice, from within the input, and from within the params section.
# If this function is called with `context == input`, we return files (if appropriate),
# so that snakemake can look for them. If `context == params` however, this is the function call
# for the params, where we return the model string if it is not a file.
def epa_ng_place_model( context ):
    if config["params"]["epa-ng"]["model-params"] == "":
        # Use the raxml best model file
        inp = "{outdir}/model/model_eval.raxml.bestModel"
        prm = ""
    elif os.path.isfile( config["params"]["epa-ng"]["model-params"] ):
        # Use the provided file
        inp = config["params"]["epa-ng"]["model-params"]
        prm = ""
    else:
        # Use the provided model string
        inp = ""
        prm = config["params"]["epa-ng"]["model-params"]

    # One of them has to be empty. This way, we can just concatenate them in the rule below.
    assert(( inp == "" ) ^ ( prm == "" ))
    if context == "input":
        return inp
    elif context == "params":
        return prm
    else:
        util.fail( "Invalid epa_ng_place_model context: " + context )
