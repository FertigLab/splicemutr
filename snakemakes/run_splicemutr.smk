import os

configfile: "config.yaml"

if not os.path.exists(config["FORMED_TRANSCRIPTS_DIR"]):
    os.mkdir(config["FORMED_TRANSCRIPTS_DIR"])
if not os.path.exists(config["COMBINE_SPLICEMUTR_OUT"]):
    os.mkdir(config["COMBINE_SPLICEMUTR_OUT"])
if not os.path.exists(config["PROCESS_PEPTIDES_OUT"]):
    os.mkdir(config["PROCESS_PEPTIDES_OUT"])
if not os.path.exists(config["MHCNUGGETS_OUT"]):
    os.mkdir(config["MHCNUGGETS_OUT"])
if not os.path.exists(config["PROCESS_BINDAFF_OUT"]):
    os.mkdir(config["PROCESS_BINDAFF_OUT"])
if not os.path.exists(config["ANALYZE_SPLICEMUTR_OUT"]):
    os.mkdir(config["ANALYZE_SPLICEMUTR_OUT"])


rule all:
    input:
        #FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds",
        #FORMED_TRANSCRIPTS_CP=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr_cp_corrected.rds",
        #OUTPUT_FILE=config["COMBINE_SPLICEMUTR_OUT"]+"/data_splicemutr_all_pep.txt",
        #OUT_FILE=config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt",
        #UNIQUE_MHC_FILE=config["GENOTYPES_DIR"]+"/class_1_HLAS.txt",
        #OUT_FILE_MHCNUGGETS=config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        #PROCESS_BINDAFF_FILES=config["PROCESS_BINDAFF_OUT"]+"/filenames.txt",
        #EXTRACT_DATA_FILE=config["PROCESS_BINDAFF_OUT"]+"/summaries.txt",
        ANALYZE_SPLICEMUTR_OUT_FILE=config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt"


rule form_transcripts:
    input:
        INTRON_FILE=config["INTRON_FILE"],
        TXDB=config["TXDB"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        SPLICEMUTR_FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"]
    output:
        FORMED_TRANSCRIPTS_DIR=config["FORMED_TRANSCRIPTS_DIR"],
        FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds"
    shell:
        """
        BSGENOME="BSgenome.Hsapiens.GENCODE.GRCh38.p10"
        OUT_PREFIX={output.FORMED_TRANSCRIPTS_DIR}/$(echo $(basename {input.INTRON_FILE}) | sed s/'.rds'/''/g)

        {input.SCRIPT_DIR}/form_transcripts.R -o $OUT_PREFIX -t {input.TXDB} -j {input.INTRON_FILE} -b $BSGENOME -f {input.SPLICEMUTR_FUNCTIONS}
        """

rule calcualte_coding_potential:
    input:
        SPLICE_FILE=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds",
        SPLICEMUTR_FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"]
    output:
        FORMED_TRANSCRIPTS_DIR=config["FORMED_TRANSCRIPTS_DIR"],
        FORMED_TRANSCRIPTS_CP=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr_cp_corrected.rds"
    shell:
        """
        TRANSCRIPT_FILE=$(echo {input.SPLICE_FILE} | sed s/'_data_splicemutr.rds'/'_sequences.fa'/g)

        {input.SCRIPT_DIR}/calc_coding_potential.R -o {output.FORMED_TRANSCRIPTS_DIR} -s {input.SPLICE_FILE} -t $TRANSCRIPT_FILE -f {input.SPLICEMUTR_FUNCTIONS}

        cd {output.FORMED_TRANSCRIPTS_DIR}
        ls $PWD/*_cp_corrected.rds > filenames_cp.txt
    
        """


rule combine_splicemutr:
    input:
        SPLICE_FILES=config["SPLICEMUTR_FILES"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"]
    output:
        OUTPUT_DIR=config["COMBINE_SPLICEMUTR_OUT"],
        OUTPUT_FILE=config["COMBINE_SPLICEMUTR_OUT"]+"/data_splicemutr_all_pep.txt"
    shell:
        """
        {input.SCRIPT_DIR}/combine_splicemutr.R -o {output.OUTPUT_DIR} -s {input.SPLICE_FILES}
        """

rule process_peptides:
    input:
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PEPTIDES=config["PROTEINS"]
    output:
        OUT_DIR=config["PROCESS_PEPTIDES_OUT"],
        OUT_FILE=config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt"
    shell:
        """
        KMER_LENGTH=9

        {input.SCRIPT_DIR}/process_peptides.py -p {input.PEPTIDES} -o {output.OUT_DIR} -k $KMER_LENGTH
        """

rule run_mhcnuggets:
    params:
        TYPE=config["TYPE"]
    input:
        INPUT_KMERS=config["INPUT_KMERS"],
        MHC_ALLELE_FILE=config["MHC_ALLELE_FILE"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"]
    output:
        OUT_DIR_MHCNUGGETS=config["MHCNUGGETS_OUT"],
        OUT_FILE_MHCNUGGETS=config["MHCNUGGETS_OUT"]+"/allele_files.txt"
    shell:
        """
        {input.SCRIPT_DIR}/runMHCnuggets.py -t {params.TYPE} -k {input.INPUT_KMERS} -m {input.MHC_ALLELE_FILE} -o {output.OUT_DIR_MHCNUGGETS}

        cd {output.OUT_DIR_MHCNUGGETS}
        ls $PWD/*_peps_9.txt > {output.OUT_FILE_MHCNUGGETS}
        """

rule process_bindaffinity:
    params:
        NUM_ALLELE_FILES=config["NUM_ALLELE_FILES"],
        KMER_LENGTH=config["KMER_LENGTH"]
    input:
        ALLELE_FILES=config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=config["PICKLE_DIR"]
    output:
        PROCESS_BINDAFF_OUT=config["PROCESS_BINDAFF_OUT"],
        PROCESS_BINDAFF_FILES=config["PROCESS_BINDAFF_OUT"]+"/filenames.txt"
    shell:
        """
        START=1
        for ((VAR=$START; VAR<={params.NUM_ALLELE_FILES}; VAR++))
        do
            ALLELE=$(sed -n ${{VAR}}p {input.ALLELE_FILES})
            BINDERS={output.PROCESS_BINDAFF_OUT}/$(echo $(basename $ALLELE) | sed 's/.txt/_filt.txt/g')
            awk -F "," '{{ if ($2 <= 500) {{ print }} }}' $ALLELE > $BINDERS

            {input.SCRIPT_DIR}/process_bindaff.py -b $BINDERS -p {input.PICKLE_DIR} -o {output.PROCESS_BINDAFF_OUT} -k {params.KMER_LENGTH}
        done

        cd {output.PROCESS_BINDAFF_OUT}
        ls $PWD/*filt.txt > filenames.txt

        """

rule extract_data:
    params:
        NUM_ALLELE_FILES=config["NUM_ALLELE_FILES"],
        KMER_SIZE_MIN=config["KMER_SIZE_MIN"],
        KMER_SIZE_MAX=config["KMER_SIZE_MAX"]
    input:
        ALLELE_FILES=config["MHC_ALLELE_FILE"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=config["PROCESS_BINDAFF_OUT"]
    output:
        EXTRACT_DATA_DIR=config["PROCESS_BINDAFF_OUT"],
        EXTRACT_DATA_FILE=config["PROCESS_BINDAFF_OUT"]+"/summaries.txt"
    shell:
        """
            START=1
            for ((VAR=$START; VAR<={params.NUM_ALLELE_FILES}; VAR++))
            do
                ALLELE=$(sed -n ${{VAR}}p {input.ALLELE_FILES})

                {input.SCRIPT_DIR}/extract_data.py -a $ALLELE -p {input.PICKLE_DIR} -b {params.KMER_SIZE_MIN} -e {params.KMER_SIZE_MAX}

            done

            cd {output.EXTRACT_DATA_DIR}
            ls $PWD/*summary.txt > summaries.txt
        """

rule analyze_splicemutr:
    params:
        SUMMARY_TYPE=config["SUMMARY_TYPE"],
        NUM_SAMPLES=config["NUM_ALLELE_FILES"]
    input:
        GENOTYPES=config["GENOTYPES_REFORMATTED"],
        SUMMARY_DIR=config["SUMMARY_DIR"],
        SPLICE_DAT_FILE=config["SPLICE_DAT_FILE"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"]
    output:
        ANALYZE_SPLICEMUTR_OUT=config["ANALYZE_SPLICEMUTR_OUT"],
        ANALYZE_SPLICEMUTR_OUT_FILE=config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt"
    shell:
        """
            START=1
            for ((VAR=$START; VAR<={params.NUM_SAMPLES}; VAR++))
            do
                $SCRIPT_DIR/analyze_splicemutr.py -g {input.GENOTYPES} -s {input.SUMMARY_DIR} -d {input.SPLICE_DAT_FILE} -o {output.OUT_DIR} -t {input.SUMMARY_TYPE} -n $VAR
            done

            cd {output.ANALYZE_SPLICEMUTR_OUT}
            ls $PWD/*_splicemutr_kmers.txt > filenames.txt
        """
