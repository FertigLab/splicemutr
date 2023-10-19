import os

configfile: "config.yaml"

os.mkdir(config["REF_DIR"])
os.chdir(config["REF_DIR"])
os.system("wget %s"%config["GTF_URL"])
os.system("wget %s"%config["FASTA_URL"])

rule get_reference_data:
    input:
        GTF_FILE_GZ=config["REF_DIR"]+"/"+config["GTF_FILE_GZ"],
        FASTA_FILE_GZ=config["REF_DIR"]+"/"+config["FASTA_FILE_GZ"]
    shell:
        """
        mkdir {input.REF_DIR}
        cd {input.REF_DIR}
        gunzip {output.GTF_FILE_GZ}
        gunzip {output.FASTA_FILE_GZ}
        """
'''
rule make_txdb:
    input:
        REF_DIR=config["REF_DIR"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
        GTF_FILE=config["GTF_FILE"]
    output:
        OUT_FILE=config["OUT_FILE"]
    shell:
        """
        conda activate miniconda3/envs/splicemutr

        {input.SPLICEMUTR_SCRIPTS}/make_txdb.R -o {input.REF_DIR}/{output.OUT_FILE} -g {input.REF_DIR}/{input.GTF_FILE}
        """
    
rule prepare_leafcutter_references:
    input:
        LEAF_DIR=config["LEAF_DIR"],
        GTF=config["GTF"]
    output:
        ANN_DIR=config["ANN_DIR"]
    shell:
        """
        {input.LEAF_DIR}/scripts/gtf_to_exons.R $GTF {ouput.ANN_DIR}/G026.exons.txt

        cd $ANN_DIR

        {input.LEAF_DIR}/leafviz/gtf2leafcutter.pl -o G026 {input.GTF)
        """

rule convert_fasta_twobit:
    input:
        REF_DIR=config["REF_DIR"],
        FA_TO_TWOBIT_EXEC=config["FA_TO_TWOBIT_EXEC"],
        FASTA_FILE=config["FASTA_FILE"]
    output:
        TWOBIT_FILE=config["TWOBIT_FILE"]
    shell:
        """
        {input.FA_TO_TWOBIT_EXEC} {input.REF_DIR}/{input.FASTA_FILE} {input.REF_DIR}/{output.TWOBIT_FILE}
        """
  
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