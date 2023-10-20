import os

configfile: "config.yaml"

if not os.path.exists(config["REF_DIR"]):
    os.mkdir(config["REF_DIR"])
os.chdir(config["REF_DIR"])
if not os.path.isfile(config["GTF_URL"]:
    os.system("wget %s"%config["GTF_URL"])
if not os.path.isfile(config["FASTA_URL"]):
    os.system("wget %s"%config["FASTA_URL"])

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

if not os.path.exists(config["ANN_DIR"]):
    os.mkdir(config["ANN_DIR"])

rule prepare_leafcutter_references:
    input:
        LEAF_DIR=config["LEAF_DIR"],
        GTF=config["GTF"]
    output:
        ANNOTATION=config["ANN_DIR"]+"/"+"G026.exons.txt"
    shell:
        """
        {input.LEAF_DIR}/scripts/gtf_to_exons.R {input.GTF} {output.ANNOTATION}

        cd $(dirname {output.ANNOTATION})

        {input.LEAF_DIR}/leafviz/gtf2leafcutter.pl -o G026 {input.GTF}
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

if os.path.exists(config["REF_DIR"]+"/"+config["BSGENOME"]):
    os.system("rm -r %s"%(config["REF_DIR"]+"/"+config["BSGENOME"]))

rule create_bsgenome:
    input:
        REF_DIR=config["REF_DIR"],
        SEED_FILE=config["SEED_FILE"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
    output:
        BSGENOME=config["REF_DIR"]+"/"+config["BSGENOME"]
    shell:
        """
        cd {input.REF_DIR}

        {input.SPLICEMUTR_SCRIPTS}/forge_bsgenome.R -s {input.SEED_FILE}
        R_BUILD_TAR=tar
        R CMD build $({output.BSGENOME})
        R CMD INSTALL {output.BSGENOME}_0.1.tar.gz
        """
