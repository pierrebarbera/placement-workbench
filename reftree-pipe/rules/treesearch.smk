# =================================================================================================
#     Helper Functions
# =================================================================================================
#
# These functions are meant to be extended via wildcards if and when more tree inference tools are added
#
def bootstrap_params( wildcards ):
    mode        = get_highest_override( ['raxml-ng', 'treesearch'], "bootstrap_metric" )
    num_trees   = get_highest_override( ['raxml-ng', 'treesearch'], "bs_trees" )
    auto_bs     = get_highest_override( ['raxml-ng', 'treesearch'], "auto_bootstrap" )

    res         = " --bs-metric {}".format(mode) if mode else ""

    if num_trees:
        if auto_bs:
            res = res + " --bs-trees autoMRE{{{}}}".format(num_trees)
        else:
            res = res + " --bs-trees {}".format(num_trees)

    return res

def model_params( wildcards, input ):
    datatype    = config["settings"]["datatype"]
    model       = get_highest_override( ['raxml-ng', 'treesearch'], "model" )

    prefix = "--model "

    if use_auto_model:
        return prefix + input.modelfile
    elif model:
        return prefix + model
    elif datatype and datatype in ['nt','aa']:
        return prefix + "GTR+G" if datatype == 'nt' else arg_string + 'LG+G'
    else:
        util.fail("'datatype' field is required, and must be either 'nt' or 'aa'.")

def model_file( wildcards ):
    if use_auto_model:
        return "{}/result/{}/{}/{}/{}/modeltest-ng/model.file".format(
            wildcards.outdir, wildcards.sample, wildcards.autoref, wildcards.aligner, wildcards.trimmer
        )
    else:
        return []

def starting_trees_params( wildcards ):
    pars_trees = get_highest_override( ['raxml-ng', 'treesearch'], "parsimony_starting_trees")
    rand_trees = get_highest_override( ['raxml-ng', 'treesearch'], "random_starting_trees")

    if pars_trees or rand_trees:
        trees = []
        if pars_trees:
            trees.append("pars{{{}}}".format(pars_trees))
        if rand_trees:
            trees.append("rand{{{}}}".format(rand_trees))
        return " --tree " + ",".join(trees)
    else:
        # nothing specified: leave it up to the tool to decide
        return ""

# =================================================================================================
#     Tree Search with RAxML-ng
# =================================================================================================

rule treesearch_raxmlng:
    group: "treesearch"
    input:
        msa = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/trimmed.afa",
        modelfile = model_file
    params:
        model           = model_params,
        starting_trees  = starting_trees_params,
        bootstrap       = bootstrap_params,
        extra           = config["params"]["raxml-ng"]["treesearch"]["extra"],
        prefix          = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/search",
        constraint      = "" if not config["data"]["constraint_tree"] else '--tree-constraint "{}"'.format(config["data"]["constraint_tree"])
    threads:
        get_threads( ['raxml-ng', 'treesearch'] )
    output:
        best_tree       = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
        best_model      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.model",
        support_tree    = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/bootstrap.newick",
        ml_trees        = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick",
        bs_trees        = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/bs_trees.newick"
    log:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/search.log"
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --all --msa {input.msa} --prefix {params.prefix}"
        " {params.model}"
        " {params.starting_trees}"
        " {params.bootstrap}"
        " --threads {threads} --redo"
        " {params.constraint}"
        " {params.extra}"
        " > {log}  2>&1"
        # symlink resulting files to be simpler to understand and conform with other methods
        " && cd $(dirname {output.best_tree})"
        " && ln -s search.raxml.bestTree best.newick"
        " && ln -s search.raxml.bestModel best.model"
        " && ln -s search.raxml.support bootstrap.newick"
        " && ln -s search.raxml.mlTrees ml_trees.newick"
        " && ln -s search.raxml.bootstraps bs_trees.newick"


# =================================================================================================
#     Consensus Tree with RAxML-ng
# =================================================================================================

rule treesearch_consensus:
    group: "treesearch"
    input:
        "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees.newick"
    output:
        mr      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/consensusTreeMR.newick",
        mre     = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/consensusTreeMRE.newick"
    params:
        prefix  = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/ml_trees"
    threads:
        get_threads( ['raxml-ng', 'treesearch'] )
    log:
        mr      = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/mr.log",
        mre     = "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/mre.log"
    conda:
        "../envs/raxml-ng.yaml"
    shell:
        "raxml-ng --consense MR  --tree {input} --prefix {params.prefix} --threads {threads} > {log.mr}  2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMR {output.mr} && "
        "raxml-ng --consense MRE --tree {input} --prefix {params.prefix} --threads {threads} > {log.mre} 2>&1 && "
        "mv {params.prefix}.raxml.consensusTreeMRE {output.mre}"

rule determine_best_run:
    """Compare all best trees of a sample and create a symlink to that folder"""
    input:
        expand(
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/tree/best.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            ),
        expand(
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMR.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            ),
        expand(  
            "{outdir}/result/{sample}/{autoref}/{aligner}/{trimmer}/raxml-ng/post/plausible.consensusTreeMRE.newick",
            autoref=autoref_list,
            aligner=aligner_list,
            trimmer=trimmer_list,
            allow_missing=True
            )
    output:
        "{outdir}/result/{sample}/best_result/raxml-ng/tree/best.newick",
        "{outdir}/result/{sample}/best_result/raxml-ng/post/plausible.consensusTreeMR.newick",
        "{outdir}/result/{sample}/best_result/raxml-ng/post/plausible.consensusTreeMRE.newick"
    script:
        "../scripts/symlink-best-result.py"

localrules: determine_best_run
