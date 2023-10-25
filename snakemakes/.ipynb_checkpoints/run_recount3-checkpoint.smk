import os

configfile: "config.yaml"

rule all:
    input:
        OUTPUT_FILE=config["JUNC_OUTPUT_DIR"]+"/filenames.txt"
        
if not os.path.exists(config["JUNC_OUTPUT_DIR"]):
    os.mkdir(config["JUNC_OUTPUT_DIR"])

rule run_recount3:
    input:
        TCGA_CANCER=config["TCGA_CANCER"],
        JUNC_OUTPUT_DIR=config["JUNC_OUTPUT_DIR"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
        SAMPLES=config["SAMPLES"]
    output:
        OUTPUT_FILE=config["JUNC_OUTPUT_DIR"]+"/filenames.txt"
    shell:
        """
        OUT_DIR={input.JUNC_OUTPUT_DIR}
        STUDY={input.TCGA_CANCER}
        SAMPLES={input.SAMPLES}
        SCRIPT_DIR={input.SPLICEMUTR_SCRIPTS}
        $SCRIPT_DIR/recount3_tcga_juncs_select.R -t $STUDY -o $OUT_DIR -s $SAMPLES

        cd {input.JUNC_OUTPUT_DIR}
        ls $PWD/*.junc > {output.OUTPUT_FILE}
        """