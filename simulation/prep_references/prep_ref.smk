import os
import platform

configfile: "config.yaml"

if not os.path.exists(config["REF_DIR"]):
    os.mkdir(config["REF_DIR"])
    os.chdir(config["REF_DIR"])
    if platform.system() == "Linux":
        os.system("wget %s"%config["GTF_URL"])
        os.system("wget %s"%config["FASTA_URL"])
    elif platform.system() == "Darwin":
        os.system("curl %s --output %s"%(config["GTF_URL"],config["GTF_FILE"]))
        os.system("curl %s --output %s"%(config["FASTA_URL"],config["FASTA_FILE_GZ"]))
if not os.path.exists(config["ANN_DIR"]):
    os.mkdir(config["ANN_DIR"])
if os.path.exists(config["BSGENOME"]):
    os.system("rm -r %s"%(config["BSGENOME"]))
if (config["REF_DIR"] in os.getcwd()):
    
rule all:
    input:
        LEAF_DIR=os.getcwd()+"/"+config["LEAF_DIR"]+"/scripts/bam2junc.sh"
        FA_TO_TWOBIT_EXEC=os.getcwd()+"/"+config["FA_TO_TWOBIT_EXEC"]
        SEED_FILE=os.getcwd()+"/"+config["SEED_FILE"]
        FASTA_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["FASTA_FILE"],
        GTF_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["GTF_FILE"],
        OUT_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["OUT_FILE"],
        ANNOTATION=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["ANN_DIR"]+"/"+"G039.exons.txt",
        TWOBIT_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["TWOBIT_FILE"],
        BSGENOME=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["BSGENOME"],

rule download_leafcutter:
    params:
        LEAF_URL=config["LEAF_URL"]
    output:
        LEAF_DIR=os.getcwd()+"/"+config["LEAF_DIR"]+"/scripts/bam2junc.sh"
    shell:
        """
        git clone {params.LEAF_URL}
        """

rule download_faToTwoBit:
    params:
        FA_TO_TWOBIT_URL=config["FA_TO_TWOBIT_URL"]
    output:
        FA_TO_TWOBIT_EXEC=os.getcwd()+"/"+config["FA_TO_TWOBIT_EXEC"]
    shell:
        """
        wget {params.FA_TO_TWOBIT_URL}
        chmod +x {output.FA_TO_TWOBIT_EXEC}
        """

rule create_description:
    output:
        SEED_FILE=os.getcwd()+"/"+config["SEED_FILE"]
    shell:
        """
        echo "Package: BSgenome.Hsapiens.GENCODE.GRCh38.p13" >> {output.SEED_FILE}
        echo "Title: Full genome sequences for Homo sapiens GRCh38.p13" >> {output.SEED_FILE}
        echo "Description: Full genome sequences for Homo sapiens (human) as provided by GENCODE (Genome Reference Consortium Human Build 38 patch release 39)" >> {output.SEED_FILE}
        echo "Version: 0.1" >> {output.SEED_FILE}
        echo "Author: Theron Palmer" >> {output.SEED_FILE}
        echo "Maintainer: Theron<tpalme15@jhmi.edu>" >> {output.SEED_FILE}
        echo "Suggests:" >> {output.SEED_FILE}
        echo "License: ?" >> {output.SEED_FILE}
        echo "organism: Homo sapiens" >> {output.SEED_FILE}
        echo "common_name: human" >> {output.SEED_FILE}
        echo "seqs_srcdir: $PWD/splicemutr_references" >> {output.SEED_FILE}
        echo "seqfile_name: GRCh38.primary_assembly.genome.2bit" >> {output.SEED_FILE}
        echo "provider: GENCODE" >> {output.SEED_FILE}
        echo "genome: GRCh38.p13" >> {output.SEED_FILE}
        echo "release_date: 12/2021" >> {output.SEED_FILE}
        echo "source_url: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/" >> {output.SEED_FILE}
        echo "organism_biocview: AnnotationData, Genetics, BSgenome, Homo_sapiens" >> {output.SEED_FILE}
        echo "BSgenomeObjname: Hsapiens" >> {output.SEED_FILE}
        """

rule get_reference_data:
    params:
        REF_DIR=os.getcwd()+"/"+config["REF_DIR"]
    input:
        FASTA_FILE_GZ=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["FASTA_FILE_GZ"]
    output:
        FASTA_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["FASTA_FILE"]
    shell:
        """
        cd {params.REF_DIR}
        gunzip -k {input.FASTA_FILE_GZ}
        """

rule make_txdb:
    input:
        REF_DIR=os.getcwd()+"/"+config["REF_DIR"],
        GTF_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["GTF_FILE"],
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    output:
        OUT_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["OUT_FILE"]
    shell:
        """

        {input.SPLICEMUTR_SCRIPTS}/make_txdb.R -o {output.OUT_FILE} -g {input.GTF_FILE}
        """

rule prepare_leafcutter_references:
    input:
        LEAF_DIR=os.getcwd()+"/"+config["LEAF_DIR"],
        GTF_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["GTF_FILE"]
    output:
        ANNOTATION=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["ANN_DIR"]+"/"+"G039.exons.txt"
    shell:
        """
        {input.LEAF_DIR}/scripts/gtf_to_exons.R {input.GTF_FILE} {output.ANNOTATION}

        cd $(dirname {output.ANNOTATION})

        {input.LEAF_DIR}/leafviz/gtf2leafcutter.pl -o G039 {input.GTF_FILE}
        """

rule convert_fasta_twobit:
    input:
        FA_TO_TWOBIT_EXEC=os.getcwd()+"/"+config["FA_TO_TWOBIT_EXEC"],
        FASTA_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["FASTA_FILE"]
    output:
        TWOBIT_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["TWOBIT_FILE"]
    shell:
        """
        {input.FA_TO_TWOBIT_EXEC} {input.FASTA_FILE} {output.TWOBIT_FILE}
        """

rule create_bsgenome:
    input:
        REF_DIR=os.getcwd()+"/"+config["REF_DIR"],
        SEED_FILE=os.getcwd()+"/"+config["SEED_FILE"],
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        TWOBIT_FILE=os.getcwd()+"/"+config["REF_DIR"]+"/"+config["TWOBIT_FILE"]
    output:
        BSGENOME=directory(os.getcwd()+"/"+config["REF_DIR"]+"/"+config["BSGENOME"])
    shell:
        """
        cd {input.REF_DIR}
        chmod +x {input.SPLICEMUTR_SCRIPTS}/*
        {input.SPLICEMUTR_SCRIPTS}/forge_bsgenome.R -s {input.SEED_FILE}
        R CMD build {output.BSGENOME}
        R CMD INSTALL {output.BSGENOME}_0.1.tar.gz
        """
