# =================================================================================================
#     Obtaining the input sequences from genbank
# =================================================================================================

rule modeltest:
	input:
        "{outdir}/result/{sample}/{aligner}/{trimmer}/trimmed.afa"
    params:
        datatype = config["params"]["pargenes"]["datatype"]
    output:
    	"{outdir}/result/{sample}/{aligner}/modeltest-ng/model.file"
    conda:
        "../envs/modeltest-ng.yaml"
    shell:
        "scripts/download_fasta.py"
