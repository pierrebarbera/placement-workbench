# =================================================================================================
#     Optional step to determine model parameters for treesearch
# =================================================================================================

rule modeltest:
    input:
       "{outdir}/result/{sample}/{aligner}/{trimmer}/trimmed.afa"
    params:
        datatype = config["settings"]["datatype"]
    threads:
        get_highest_override( "modeltestng", "threads" )
    output:
        "{outdir}/result/{sample}/{aligner}/{trimmer}/modeltest-ng/model.file"
    log:
        "{outdir}/result/{sample}/{aligner}/{trimmer}/modeltest-ng/modeltest.log"
    conda:
        "../envs/modeltest-ng.yaml"
    script:
        "../scripts/modeltest-ng.py"
