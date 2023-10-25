import os

configfile: "config.yaml"

if not os.path.exists(config["FORMED_TRANSCRIPTS_DIR"]):
    os.mkdir(config["FORMED_TRANSCRIPTS_DIR"])

rule all:
    input:
        
    
rule form_transcripts:
    input:
        INTRON_FILE=config["INTRON_FILE"],
        TXDB=config["TXDB"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        SPLICEMUTR_FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"]
    output:
        FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds"
    shell:
        """
        BSGENOME="BSgenome.Hsapiens.GENCODE.GRCh38.p10"
        OUT_PREFIX={output.FORMED_TRANSCRIPTS_DIR}/$(echo $(basename {input.INTRON_FILE}) | sed s/'.rds'/''/g)

        {input.SCRIPT_DIR}/form_transcripts.R -o $OUT_PREFIX -t {input.TXDB} -j {input.INTRON_FILE} -b $BSGENOME
        """