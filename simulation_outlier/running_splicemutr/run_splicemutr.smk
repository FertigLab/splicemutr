import os

configfile: "config.yaml"

if not os.path.exists(os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]):
    os.mkdir(os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"])
if not os.path.exists(os.getcwd()+"/"+config["COMBINE_SPLICEMUTR_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["COMBINE_SPLICEMUTR_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["PROCESS_PEPTIDES_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["PROCESS_PEPTIDES_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["MHCNUGGETS_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["MHCNUGGETS_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"])
if not os.path.exists(os.getcwd()+"/"+config["CREATE_JUNC_EXPRESSION_OUT"]):
    os.mkdir(os.getcwd()+"/"+config["CREATE_JUNC_EXPRESSION_OUT"])

rule all:
    input:
        CLASS_1_ALLELES=os.getcwd()+"/class1_alleles.txt",
        FORMED_TRANSCRIPTS=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/comparison_juncs_linear_data_splicemutr.rds",
        FORMED_TRANSCRIPTS_CP=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/filenames_cp.txt",
        OUTPUT_FILE=os.getcwd()+"/"+config["COMBINE_SPLICEMUTR_OUT"]+"/proteins.txt",
        SPLICE_DAT_FILE=os.getcwd()+"/"+config["SPLICE_DAT_FILE"],
        OUT_FILE=os.getcwd()+"/"+config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt",
        OUT_FILE_MHCNUGGETS=os.getcwd()+"/"+config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        PROCESS_BINDAFF_FILES=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/filenames.txt",
        EXTRACT_DATA_FILE=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/summaries.txt",
        PROTEIN_FASTA=os.getcwd()+"/"+config["PROTEIN_FASTA"],
        ANALYZE_SPLICEMUTR_OUT_FILE=os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt",
        COMPARISONS_FILE=os.getcwd()+"/"+config["COMPARISONS_FILE"],
        COMPARISONS_OUT_FILE=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/create_comparisons_out_cp/splice_dat_data.rds",
        KMER_COUNTS_FILE=os.getcwd()+"/"+config["KMER_COUNTS_FILE"],
        CREATE_JUNC_EXPRESSION_FILE=config["CREATE_JUNC_EXPRESSION_OUT"]+"/junc_expr_combined_vst.rds",
        SPLICING_ANTIGENICITY_FILE=os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"]+"/filenames.txt"

rule create_class_1_alleles_file:
    output:
        CLASS_1_ALLELES=os.getcwd()+"/class1_alleles.txt"
    shell:
        """
        echo "HLA-A01-01" >> {output.CLASS_1_ALLELES}
        """

rule form_transcripts:
    input:
        INTRON_FILE=os.getcwd()+"/"+config["INTRON_FILE"],
        TXDB=os.getcwd()+"/"+config["TXDB"],
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        SPLICEMUTR_FUNCTIONS=os.getcwd()+"/"+config["SPLICEMUTR_FUNCTIONS"],
        FORMED_TRANSCRIPTS_DIR=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]
    output:
        FORMED_TRANSCRIPTS=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/comparison_juncs_linear_data_splicemutr.rds"
    shell:
        """
        BSGENOME="BSgenome.Hsapiens.GENCODE.GRCh38.p13"
        OUT_PREFIX={input.FORMED_TRANSCRIPTS_DIR}/$(echo $(basename {input.INTRON_FILE}) | sed s/'.rds'/''/g)

        {input.SCRIPT_DIR}/form_transcripts.R -o $OUT_PREFIX -t {input.TXDB} -j {input.INTRON_FILE} -b $BSGENOME -f {input.SPLICEMUTR_FUNCTIONS}
        """

rule calcualte_coding_potential:
    input:
        SPLICE_FILE=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/comparison_juncs_linear_data_splicemutr.rds",
        SPLICEMUTR_FUNCTIONS=os.getcwd()+"/"+config["SPLICEMUTR_FUNCTIONS"],
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        FORMED_TRANSCRIPTS_DIR=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]
    output:
        SPLICE_FILES=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/filenames_cp.txt"
    shell:
        """
        TRANSCRIPT_FILE=$(echo {input.SPLICE_FILE} | sed s/'_data_splicemutr.rds'/'_sequences.fa'/g)

        {input.SCRIPT_DIR}/calc_coding_potential.R -o {input.FORMED_TRANSCRIPTS_DIR} -s {input.SPLICE_FILE} -t $TRANSCRIPT_FILE -f {input.SPLICEMUTR_FUNCTIONS}

        cd {input.FORMED_TRANSCRIPTS_DIR}
        ls $PWD/*_cp_corrected.rds > filenames_cp.txt
    
        """

rule combine_splicemutr:
    input:
        SPLICE_FILES=os.getcwd()+"/"+config["FORMED_TRANSCRIPTS_DIR"]+"/filenames_cp.txt",
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        OUTPUT_DIR=os.getcwd()+"/"+config["COMBINE_SPLICEMUTR_OUT"]
    output:
        OUTPUT_FILE=os.getcwd()+"/"+config["COMBINE_SPLICEMUTR_OUT"]+"/proteins.txt",
        SPLICE_DAT_FILE=os.getcwd()+"/"+config["SPLICE_DAT_FILE"]
    shell:
        """
        {input.SCRIPT_DIR}/combine_splicemutr.R -o {input.OUTPUT_DIR} -s {input.SPLICE_FILES}
        """

rule process_peptides:
    input:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_PYTHON"],
        PEPTIDES=os.getcwd()+"/"+config["PROTEINS"],
        OUT_DIR=os.getcwd()+"/"+config["PROCESS_PEPTIDES_OUT"]
    output:
        OUT_FILE=os.getcwd()+"/"+config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt"
    shell:
        """
        KMER_LENGTH=9

        {input.SCRIPT_DIR}/process_peptides.py -p {input.PEPTIDES} -o {input.OUT_DIR} -k $KMER_LENGTH
        """

rule run_mhcnuggets:
    params:
        TYPE=config["TYPE"]
    input:
        INPUT_KMERS=os.getcwd()+"/"+config["INPUT_KMERS"],
        CLASS_1_ALLELES=os.getcwd()+"/class1_alleles.txt",
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_PYTHON"],
        OUT_DIR_MHCNUGGETS=os.getcwd()+"/"+config["MHCNUGGETS_OUT"]
    output:
        OUT_FILE_MHCNUGGETS=os.getcwd()+"/"+config["MHCNUGGETS_OUT"]+"/allele_files.txt"
    shell:
        """
        {input.SCRIPT_DIR}/runMHCnuggets.py -t {params.TYPE} -k {input.INPUT_KMERS} -m {input.CLASS_1_ALLELES} -o {input.OUT_DIR_MHCNUGGETS}

        cd {input.OUT_DIR_MHCNUGGETS}
        ls $PWD/*_peps_9.txt > {output.OUT_FILE_MHCNUGGETS}
        """

rule process_bindaffinity:
    params:
        NUM_ALLELE_FILES=config["NUM_ALLELE_FILES"],
        KMER_LENGTH=config["KMER_LENGTH"]
    input:
        ALLELE_FILES=os.getcwd()+"/"+config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=os.getcwd()+"/"+config["PICKLE_DIR"],
        PROCESS_BINDAFF_OUT=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]
    output:
        PROCESS_BINDAFF_FILES=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/filenames.txt"
    shell:
        """
        START=1
        for ((VAR=$START; VAR<={params.NUM_ALLELE_FILES}; VAR++))
        do
            ALLELE=$(sed -n ${{VAR}}p {input.ALLELE_FILES})
            BINDERS={input.PROCESS_BINDAFF_OUT}/$(echo $(basename $ALLELE) | sed 's/.txt/_filt.txt/g')
            awk -F "," '{{ if ($2 <= 500) {{ print }} }}' $ALLELE > $BINDERS

            {input.SCRIPT_DIR}/process_bindaff.py -b $BINDERS -p {input.PICKLE_DIR} -o {input.PROCESS_BINDAFF_OUT} -k {params.KMER_LENGTH}
        done

        cd {input.PROCESS_BINDAFF_OUT}
        ls $PWD/*filt.txt > filenames.txt

        """

rule extract_data:
    params:
        NUM_ALLELE_FILES=config["NUM_ALLELE_FILES"],
        KMER_SIZE_MIN=config["KMER_SIZE_MIN"],
        KMER_SIZE_MAX=config["KMER_SIZE_MAX"]
    input:
        CLASS_1_ALLELES=os.getcwd()+"/class1_alleles.txt",
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"],
        EXTRACT_DATA_DIR=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"],
        PROCESS_BINDAFF_FILES=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/filenames.txt"
    output:
        EXTRACT_DATA_FILE=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/summaries.txt"
    shell:
        """
            START=1
            for ((VAR=$START; VAR<={params.NUM_ALLELE_FILES}; VAR++))
            do
                ALLELE=$(sed -n ${{VAR}}p {input.CLASS_1_ALLELES})

                {input.SCRIPT_DIR}/extract_data.py -a $ALLELE -p {input.PICKLE_DIR} -b {params.KMER_SIZE_MIN} -e {params.KMER_SIZE_MAX}

            done

            cd {input.EXTRACT_DATA_DIR}
            ls $PWD/*summary.txt > summaries.txt
        """

rule obtain_reference_fasta:
    params:
        PROTEIN_FASTA_HTML=config["PROTEIN_FASTA_HTML"]
    output:
        PROTEIN_FASTA=os.getcwd()+"/"+config["PROTEIN_FASTA"]
    shell:
        """
        wget {params.PROTEIN_FASTA_HTML}
        gunzip *fa.gz
        """

rule create_genotypes:
    input:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    output:
        GENOTYPES=os.getcwd()+"/genotypes.rds"
    shell:
        """
            chmod +x {input.SCRIPT_DIR}/*
            {input.SCRIPT_DIR}/create_sim_genotypes.R
        """

rule analyze_splicemutr:
    params:
        SUMMARY_DIR=os.getcwd()+"/"+config["SUMMARY_DIR"],
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        OUT_DIR=os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"]
    input:
        OUT_FILE_MHCNUGGETS=os.getcwd()+"/"+config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        PROCESS_BINDAFF_FILES=os.getcwd()+"/"+config["PROCESS_BINDAFF_OUT"]+"/filenames.txt",
        GENOTYPES=os.getcwd()+"/"+config["GENOTYPES"],
        SPLICE_DAT_FILE=os.getcwd()+"/"+config["SPLICE_DAT_FILE"],
        COUNTS_FILE=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/sample_01.filt.junc",
        ANALYZE_SPLICEMUTR_OUT=os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"],
        PROTEIN_FASTA=os.getcwd()+"/"+config["PROTEIN_FASTA"],
    output:
        ANALYZE_SPLICEMUTR_OUT_FILE=os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt"
    shell:
        """
            {params.SCRIPT_DIR}/valsamo_analyze_splicemutr.R -g {input.GENOTYPES} -s {params.SUMMARY_DIR} -d {input.SPLICE_DAT_FILE} -c {input.COUNTS_FILE} -o {params.OUT_DIR} -r {input.PROTEIN_FASTA}
        """

rule create_comparisons_file:
    params:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    output:
        COMPARISONS_FILE=os.getcwd()+"/"+config["COMPARISONS_FILE"]
    shell:
        """
            chmod +x {params.SCRIPT_DIR}/*
            {params.SCRIPT_DIR}/create_comparisons_file.R
        """

rule create_comparisons:
    params:
        COMP_NUM=1,
        OUT_DIR=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/create_comparisons_out_cp",
        LEAF_DIR=config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT",
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    input:
        COMPS_JUNCS_FILE=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/comparison_juncs.rds",
        SPLICE_DAT_FILE=os.getcwd()+"/"+config["SPLICE_DAT_FILE"],
        COMPARISONS_FILE=os.getcwd()+"/"+config["COMPARISONS_FILE"],
        JUNC_DIR=os.getcwd()+"/"+config["ANALYZE_SPLICEMUTR_OUT"],
    output:
        COMPARISONS_OUT_FILE=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/create_comparisons_out_cp/splice_dat_data.rds",
        KMER_COUNTS_FILE=os.getcwd()+"/"+config["KMER_COUNTS_FILE"]
    shell:
        """
            mkdir -p {params.OUT_DIR}
            {input.SCRIPT_DIR}/create_outlier_single_comparisons.R -c {input.COMPS_JUNCS_FILE} -d {input.SPLICE_DAT_FILE} -e {input.COMPARISONS_FILE} -n {params.COMP_NUM} -j {params.JUNC_DIR} -o {params.OUT_DIR} -l {params.LEAF_DIR}
        """

rule create_junction_expression:
    input:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        JUNC_FILES=os.getcwd()+"/"+config["JUNCFILES"],
        CREATE_JUNC_EXPRESSION_OUT=os.getcwd()+"/"+config["CREATE_JUNC_EXPRESSION_OUT"]
    output:
        CREATE_JUNC_EXPRESSION_FILE=config["CREATE_JUNC_EXPRESSION_OUT"]+"/junc_expr_combined_vst.rds"
    shell:
        """
        {input.SCRIPT_DIR}/create_junc_expr.R -f {input.JUNC_FILES} -o {input.CREATE_JUNC_EXPRESSION_OUT}
        """

rule calculate_gene_metric:
    params:
        CREATE_SPLICING_ANTIGENICITY_OUT=os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"]
    input:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        SPLICE_DAT_FILE=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/create_comparisons_out_cp/splice_dat_data.rds",
        KMER_COUNTS_FILE=os.getcwd()+"/"+config["KMER_COUNTS_FILE"],
        CREATE_JUNC_EXPRESSION_FILE=config["CREATE_JUNC_EXPRESSION_OUT"]+"/junc_expr_combined_vst.rds",
        CREATE_SPLICING_ANTIGENICITY_OUT=os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"]
    output:
        SPLICING_ANTIGENICITY_FILE=os.getcwd()+"/"+config["CREATE_SPLICING_ANTIGENICITY_OUT"]+"/filenames.txt"
    shell:
        """
        mkdir -p {input.CREATE_SPLICING_ANTIGENICITY_OUT}
        OUT_PREFIX=/dcs04/fertig/data/theron/share/run_20230320/calc_gene_metric_out/baseline_POST

        {input.SCRIPT_DIR}/calc_gene_metric_len_norm.R -s {input.SPLICE_DAT_FILE} -k {input.KMER_COUNTS} -j {input.CREATE_JUNC_EXPRESSION_FILE} -o {params.CREATE_SPLICING_ANTIGENICITY_OUT}/splicemutr

        cd {input.CREATE_SPLICING_ANTIGENICITY_OUT}
        ls $PWD/*_gene_metric_mean_len_norm_no_gene_norm_tumor.rds > filenames.txt
        """
