import os

configfile: "config.yaml"

if not os.path.exists(config["REF_DIR"]):
    os.mkdir(config["REF_DIR"])
    os.chdir(config["REF_DIR"])
    os.system("wget %s"%config["GTF_URL"])
    os.system("wget %s"%config["FASTA_URL"])
if not os.path.exists(config["ANN_DIR"]):
    os.mkdir(config["ANN_DIR"])
if os.path.exists(config["REF_DIR"]+"/"+config["BSGENOME"]):
    os.system("rm -r %s"%(config["REF_DIR"]+"/"+config["BSGENOME"]))

rule all:
    input:
        FASTA_FILE=config["REF_DIR"]+"/"+config["FASTA_FILE"],
        GTF_FILE=config["REF_DIR"]+"/"+config["GTF_FILE"],
        OUT_FILE=config["REF_DIR"]+"/"+config["OUT_FILE"],
        ANNOTATION=config["ANN_DIR"]+"/"+"G026.exons.txt",
        TWOBIT_FILE=config["REF_DIR"]+"/"+config["TWOBIT_FILE"],
        BSGENOME=config["REF_DIR"]+"/"+config["BSGENOME"]

rule get_reference_data:
    input:
        REF_DIR=config["REF_DIR"],
        FASTA_FILE_GZ=config["REF_DIR"]+"/"+config["FASTA_FILE_GZ"]
    output:
        FASTA_FILE=config["REF_DIR"]+"/"+config["FASTA_FILE"],
    shell:
        """
        cd {input.REF_DIR}
        gunzip -k {input.FASTA_FILE_GZ}
        """

rule make_txdb:
    input:
        REF_DIR=config["REF_DIR"],
        GTF_FILE=config["REF_DIR"]+"/"+config["GTF_FILE"]
    output:
        OUT_FILE=config["REF_DIR"]+"/"+config["OUT_FILE"]
    shell:
        """

        {input.SPLICEMUTR_SCRIPTS}/make_txdb.R -o {output.OUT_FILE} -g {input.GTF_FILE}
        """

rule prepare_leafcutter_references:
    input:
        LEAF_DIR=config["LEAF_DIR"],
        GTF_FILE=config["REF_DIR"]+"/"+config["GTF_FILE"]
    output:
        ANNOTATION=config["ANN_DIR"]+"/"+"G026.exons.txt"
    shell:
        """
        {input.LEAF_DIR}/scripts/gtf_to_exons.R {input.GTF_FILE} {output.ANNOTATION}

        cd $(dirname {output.ANNOTATION})

        {input.LEAF_DIR}/leafviz/gtf2leafcutter.pl -o G026 {input.GTF_FILE}
        """

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

rule create_bsgenome:
    input:
        REF_DIR=config["REF_DIR"],
        SEED_FILE=config["SEED_FILE"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"]
    output:
        BSGENOME=config["REF_DIR"]+"/"+config["BSGENOME"]
    shell:
        """
        cd {input.REF_DIR}

        {input.SPLICEMUTR_SCRIPTS}/forge_bsgenome.R -s {input.SEED_FILE}
        R CMD build {output.BSGENOME}
        R CMD INSTALL {output.BSGENOME}_0.1.tar.gz
        """
