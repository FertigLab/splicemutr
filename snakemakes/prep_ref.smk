import os

configfile: "config.yaml"

if not os.path.exists(config["REF_DIR"]):
    os.mkdir(config["REF_DIR"])
os.chdir(config["REF_DIR"])
if not os.path.isfile(config["GTF_URL"]:
    os.system("wget %s"%config["GTF_URL"])
if not os.path.isfile(config["FASTA_URL"])
    os.system("wget %s"%config["FASTA_URL"])
    
'''
rule get_reference_data:
    input:
        REF_DIR=config["REF_DIR"]
        GTF_FILE_GZ=config["REF_DIR"]+"/"+config["GTF_FILE_GZ"],
        FASTA_FILE_GZ=config["REF_DIR"]+"/"+config["FASTA_FILE_GZ"]
    shell:
        """
        cd {input.REF_DIR}
        gunzip {input.GTF_FILE_GZ}
        gunzip {input.FASTA_FILE_GZ}
        """

rule make_txdb:
    input:
        REF_DIR=config["REF_DIR"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
        GTF_FILE=config["REF_DIR"]+"/"+config["GTF_FILE"]
    output:
        OUT_FILE=config["REF_DIR"]+"/"+config["OUT_FILE"]
    shell:
        """

        {input.SPLICEMUTR_SCRIPTS}/make_txdb.R -o {output.OUT_FILE} -g {input.GTF_FILE}
        """
'''

rule convert_fasta_twobit:
    input:
        FA_TO_TWOBIT_EXEC=config["FA_TO_TWOBIT_EXEC"],
        FASTA_FILE=config["REF_DIR"]+"/"+config["FASTA_FILE"]
    output:
        TWOBIT_FILE=config["REF_DIR"]+"/"+config["TWOBIT_FILE"]
    shell:
        """
        {input.FA_TO_TWOBIT_EXEC} {input.FASTA_FILE} {output.TWOBIT_FILE}
        """
'''

rule create_bsgenome:
    input:
        REF_DIR=config["REF_DIR"],
        SEED_FILE=config["SEED_FILE"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
        BSGENOME=config["BSGENOME"]
    shell:
        """
        conda activate splicemutr
    
        cd {input.REF_DIR}
    
        {input.SPLICEMUTR_SCRIPTS}/forge_bsgenome.R -s {input.SEED_FILE}
    
        R CMD BUILD {input.BSGENOME}
        R CMD check {input.BSGENOME}.tar.gz
        R CMD INSTALL {input.BSGENOME}.tar.gz
        """
'''