# Rule to run gappa phat algorithm from a raw database of potential reference sequences
rule gappa_phat_raw:
    input:
        # the taxonomy file is expected to be passed via the config, for now.
        taxonomy_file=get_taxonomy_file,
        # special call to get_fasta, as we are calling it from the optional PhAT step
        sequence_file=get_fasta( wildcards, True )
    output:
        "{outdir}/result/{sample}/phat/ref_candidates.fa"
    params:
        # the target size is a required input for the PhAT algorithm
        target_size = config["settings"]["target_taxa_number"],
        extra       = config["params"]["gappa"]["phat"]["extra"]
    log:
        "{outdir}/result/{sample}/phat/{sample}.log"
    conda:
        "../envs/gappa.yaml"
    script:
        "scripts/gappa-phat.py"
