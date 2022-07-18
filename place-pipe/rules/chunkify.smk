# =================================================================================================
#     Chunkify and Unchunkify
# =================================================================================================

# We want the sample names of the output files to fit the naming given by the user.
# In case our input for the samples is a table with two columns, the sample names (first column)
# can be different from the file names. However, the gappa chunkify command uses file names as
# sample names. For now, it is easiest to just symlink to the files to get a list of properly
# named files. In the future, we might add an option to rename samples in gappa chunkify,
# to make this step a bit less convoluted...
rule chunkify_sample_prep:
    input:
        fasta = get_sample_fasta
    params:
        rel_input = relative_input_path,
    output:
        "{outdir}/chunkify/samples/{sample}.fasta"
    shell:
        "ln -s {params.rel_input} {output}"

# No need to execute this on the cluster computed nodes.
localrules: chunkify_sample_prep

# The rule to chunkify input fasta samples (query sequences) into chunks of equal size without
# duplicate sequences.
# We need a snakemake checkpoint here, because we cannot predict the number of chunks being produced.
checkpoint chunkify:
    input:
        # Request renamed samples, using the rule above, to get chunkify to use proper sample names.
        expand( "{outdir}/chunkify/samples/{sample}.fasta", outdir=outdir, sample=sample_names )
    output:
        abundances  = expand(   "{outdir}/chunkify/abundances/abundances_{sample}.json",
                                outdir=outdir,
                                sample=sample_names
                            ),
        chunks_dir  = directory( f"{outdir}/chunkify/chunks" )
    params:
        chunks_dir      = lambda wildcards: f"{outdir}/chunkify/chunks",
        abundances_dir  = lambda wildcards: f"{outdir}/chunkify/abundances",
        hashfunction    = config["params"]["chunkify"]["hash-function"],
        minabun         = config["params"]["chunkify"]["min-abundance"],
        chunksize       = config["params"]["chunkify"]["chunk-size"]
    log:
        expand("{outdir}/chunkify/chunkify_{sample}.log", outdir=outdir, sample=sample_names)
    conda:
        "../envs/gappa.yaml"
    shell:
        "mkdir -p {params.chunks_dir} ;"
        " mkdir -p {params.abundances_dir} ;"
        " gappa prepare chunkify"
        " --fasta-path {input}"
        " --chunks-out-dir {params.chunks_dir}"
        " --abundances-out-dir {params.abundances_dir}"
        " --hash-function {params.hashfunction}"
        " --min-abundance {params.minabun}"
        " --chunk-size {params.chunksize}"
        " > {log} 2>&1"

# Following the documentation tutorial here:
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution
def aggregate_chunkify_chunks(wildcards):
    # Wildcards are ignored here, as the chunkify process is for the whole dataset,
    # so we do not have any wildcards to take into account;
    # but still have to use them in the argument list,
    # because otherwise snakemake complains.
    chunks     = checkpoints.chunkify.get().output["chunks_dir"]
    # abundances = checkpoints.chunkify.get().output[1]
    return expand(
        "{outdir}/chunkify/placed/{chunk}.jplace",
        chunk = glob_wildcards( os.path.join(chunks, "{chunk}.fasta")).chunk,
        outdir=outdir
    )

rule unchunkify:
    input:
        aggregate_chunkify_chunks,
        expand("{outdir}/chunkify/abundances/abundances_{sample}.json", outdir=outdir, sample=sample_names)
    output:
        protected( expand( "{outdir}/placed/{sample}.jplace", outdir=outdir, sample=sample_names ))
    params:
        hash_function = config["params"]["chunkify"]["hash-function"]
    log:
        expand("{outdir}/chunkify/unchunkify_{sample}.log", outdir=outdir, sample=sample_names)
    conda:
        "../envs/gappa.yaml"
    shell:
        "gappa prepare unchunkify"
        " --abundances-path {outdir}/chunkify/abundances"
        " --chunk-file-expression {outdir}/chunkify/placed/chunk_@.jplace"
        " --hash-function {params.hash_function}"
        " --out-dir {outdir}/placed"
        " > {log} 2>&1"
