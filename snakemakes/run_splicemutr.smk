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
if not os.path.exists(config["CREATE_SPLICING_ANTIGENICITY_OUT"]):
    os.mkdir(config["CREATE_SPLICING_ANTIGENICITY_OUT"])
if not os.path.exists(config["JUNC_OUTPUT_DIR"]):
    os.mkdir(config["JUNC_OUTPUT_DIR"])
if not os.path.exists(config["JUNC_DIR"]):
    os.mkdir(config["JUNC_DIR"])
if not os.path.exists(config["CREATE_JUNC_EXPRESSION_OUT"]):
    os.mkdir(config["CREATE_JUNC_EXPRESSION_OUT"])
if not os.path.exists(config["KMER_COUNTS_OUT"]):
    os.mkdir(config["KMER_COUNTS_OUT"])

'''
rule all:
    input:
        OUTPUT_FILE_RECOUNT3=config["JUNC_OUTPUT_DIR"]+"/filenames.txt",
        RDATA=config["JUNC_DIR"]+"/data.Rdata",
        OUT_FILE_FINAL=config["JUNC_DIR"]+"/CHOL_introns.rds",
        FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds",
        FORMED_TRANSCRIPTS_CP=config["FORMED_TRANSCRIPTS_DIR"]+"/filenames_cp.txt",
        OUTPUT_FILE=config["COMBINE_SPLICEMUTR_OUT"]+"/proteins.txt",
        SPLICE_DAT_FILE=config["SPLICE_DAT_FILE"],
        OUT_FILE=config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt",
        OUT_FILE_MHCNUGGETS=config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        PROCESS_BINDAFF_FILES=config["PROCESS_BINDAFF_OUT"]+"/filenames.txt",
        EXTRACT_DATA_FILE=config["PROCESS_BINDAFF_OUT"]+"/summaries.txt",
        ANALYZE_SPLICEMUTR_OUT_FILE=config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt",
        KMER_COUNTS_FILE=config["KMER_COUNTS_OUT"]+"/all_kmers_counts.txt",
        CREATE_JUNC_EXPRESSION_FILE=config["CREATE_JUNC_EXPRESSION_OUT"]+"/junc_expr_combined_vst.rds",
        SPLICING_ANTIGENICITY_FILE=config["CREATE_SPLICING_ANTIGENICITY_OUT"]+"/filenames.txt"
'''

rule run_recount3:
    input:
        JUNC_OUTPUT_DIR=config["JUNC_OUTPUT_DIR"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"]
    output:
        OUTPUT_FILE_RECOUNT3=config["JUNC_OUTPUT_DIR"]+"/filenames.txt",
    shell:
        """
        OUT_DIR={input.JUNC_OUTPUT_DIR}
        STUDY="CHOL"
        SAMPLES="all"
        SCRIPT_DIR={input.SPLICEMUTR_SCRIPTS}
        $SCRIPT_DIR/recount3_tcga_juncs_select.R -t $STUDY -o $OUT_DIR -s $SAMPLES

        cd {input.JUNC_OUTPUT_DIR}
        ls $PWD/*.junc > {output.OUTPUT_FILE_RECOUNT3}
        """

rule running_leafcutter:
    input:
        JUNC_DIR=config["JUNC_DIR"],
        JUNCFILE_FILENAMES=config["JUNCFILE_FILENAMES"],
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"],
        LEAFCUTTER_SCRIPTS=config["LEAFCUTTER_SCRIPTS"],
        REF_DIR=config["ANN_DIR"],
        LEAFVIZ_DIR=config["LEAFVIZ_DIR"],
        GROUPS_FILE=config["GROUPS_FILE"],
        LEAFCUTTER_PYTHON=config["LEAFCUTTER_PYTHON"],
    output:
        RDATA=config["JUNC_DIR"]+"/data.Rdata"
    shell:
        """
        echo "leafcutter_cluster_regtools"
        python {input.LEAFCUTTER_PYTHON}/splicemutr_leafcutter_cluster_regtools.py -j {input.JUNCFILE_FILENAMES} -r {input.JUNC_DIR} -o data -l 500000

        echo "leafcutter_ds"
        {input.LEAFCUTTER_SCRIPTS}/leafcutter_ds.R --num_threads 1 --exon_file={input.REF_DIR}/G026.exons.txt.gz -o {input.JUNC_DIR}/leafcutter_ds {input.JUNC_DIR}/data_perind_numers.counts.gz {input.GROUPS_FILE}

        echo "prepare_results"
        {input.LEAFVIZ_DIR}/prepare_results.R -o {output.RDATA} -m {input.GROUPS_FILE} {input.JUNC_DIR}/data_perind_numers.counts.gz {input.JUNC_DIR}/leafcutter_ds_cluster_significance.txt {input.JUNC_DIR}/leafcutter_ds_effect_sizes.txt {input.REF_DIR}/G026
        """

rule save_introns:
    input:
        RDATA=config["JUNC_DIR"]+"/data.Rdata",
        SPLICEMUTR_SCRIPTS=config["SPLICEMUTR_SCRIPTS"]  
    output:
        OUT_FILE_FINAL=config["JUNC_DIR"]+"/CHOL_introns.rds"
    shell:
        """
        {input.SPLICEMUTR_SCRIPTS}/save_introns.R -i {input.RDATA} -o CHOL
        """

rule form_transcripts:
    input:
        INTRON_FILE=config["INTRON_FILE"],
        TXDB=config["TXDB"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        SPLICEMUTR_FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"],
        FORMED_TRANSCRIPTS_DIR=config["FORMED_TRANSCRIPTS_DIR"]
    output:
        FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds"
    shell:
        """
        BSGENOME="BSgenome.Hsapiens.GENCODE.GRCh38.p10"
        OUT_PREFIX={input.FORMED_TRANSCRIPTS_DIR}/$(echo $(basename {input.INTRON_FILE}) | sed s/'.rds'/''/g)

        {input.SCRIPT_DIR}/form_transcripts.R -o $OUT_PREFIX -t {input.TXDB} -j {input.INTRON_FILE} -b $BSGENOME -f {input.SPLICEMUTR_FUNCTIONS}
        """

rule calcualte_coding_potential:
    input:
        SPLICE_FILE=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds",
        SPLICEMUTR_FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        FORMED_TRANSCRIPTS_DIR=config["FORMED_TRANSCRIPTS_DIR"]
    output:
        FORMED_TRANSCRIPTS_CP=config["FORMED_TRANSCRIPTS_DIR"]+"/filenames_cp.txt"
    shell:
        """
        TRANSCRIPT_FILE=$(echo {input.SPLICE_FILE} | sed s/'_data_splicemutr.rds'/'_sequences.fa'/g)

        {input.SCRIPT_DIR}/calc_coding_potential.R -o {input.FORMED_TRANSCRIPTS_DIR} -s {input.SPLICE_FILE} -t $TRANSCRIPT_FILE -f {input.SPLICEMUTR_FUNCTIONS}

        cd {input.FORMED_TRANSCRIPTS_DIR}
        ls $PWD/*_cp_corrected.rds > filenames_cp.txt
    
        """

rule combine_splicemutr:
    input:
        SPLICE_FILES=config["SPLICEMUTR_FILES"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        OUTPUT_DIR=config["COMBINE_SPLICEMUTR_OUT"]
    output:
        OUTPUT_FILE=config["COMBINE_SPLICEMUTR_OUT"]+"/proteins.txt",
        SPLICE_DAT_FILE=config["SPLICE_DAT_FILE"]
    shell:
        """
        {input.SCRIPT_DIR}/combine_splicemutr.R -o {input.OUTPUT_DIR} -s {input.SPLICE_FILES}
        """

rule process_peptides:
    input:
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PEPTIDES=config["PROTEINS"],
        OUT_DIR=config["PROCESS_PEPTIDES_OUT"]
    output:
        OUT_FILE=config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt"
    shell:
        """
        KMER_LENGTH=9

        {input.SCRIPT_DIR}/process_peptides.py -p {input.PEPTIDES} -o {input.OUT_DIR} -k $KMER_LENGTH
        """

rule run_mhcnuggets:
    params:
        TYPE=config["TYPE"]
    input:
        INPUT_KMERS=config["INPUT_KMERS"],
        MHC_ALLELE_FILE=config["MHC_ALLELE_FILE"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        OUT_DIR_MHCNUGGETS=config["MHCNUGGETS_OUT"]
    output:
        OUT_FILE_MHCNUGGETS=config["MHCNUGGETS_OUT"]+"/allele_files.txt"
    shell:
        """
        {input.SCRIPT_DIR}/runMHCnuggets.py -t {params.TYPE} -k {input.INPUT_KMERS} -m {input.MHC_ALLELE_FILE} -o {input.OUT_DIR_MHCNUGGETS}

        cd {input.OUT_DIR_MHCNUGGETS}
        ls $PWD/*_peps_9.txt > {output.OUT_FILE_MHCNUGGETS}
        """

rule process_bindaffinity:
    params:
        NUM_ALLELE_FILES=config["NUM_ALLELE_FILES"],
        KMER_LENGTH=config["KMER_LENGTH"]
    input:
        ALLELE_FILES=config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=config["PICKLE_DIR"],
        PROCESS_BINDAFF_OUT=config["PROCESS_BINDAFF_OUT"]
    output:
        PROCESS_BINDAFF_FILES=config["PROCESS_BINDAFF_OUT"]+"/filenames.txt"
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
        KMER_SIZE_MAX=config["KMER_SIZE_MAX"],
    input:
        ALLELE_FILES=config["MHC_ALLELE_FILE"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        PICKLE_DIR=config["PROCESS_BINDAFF_OUT"],
        EXTRACT_DATA_DIR=config["PROCESS_BINDAFF_OUT"]
    output:
        EXTRACT_DATA_FILE=config["PROCESS_BINDAFF_OUT"]+"/summaries.txt"
    shell:
        """
            START=1
            for ((VAR=$START; VAR<={params.NUM_ALLELE_FILES}; VAR++))
            do
                ALLELE=$(sed -n ${{VAR}}p {input.ALLELE_FILES})

                {input.SCRIPT_DIR}/extract_data.py -a $ALLELE -p {input.PICKLE_DIR} -b {params.KMER_SIZE_MIN} -e {params.KMER_SIZE_MAX}

            done

            cd {input.EXTRACT_DATA_DIR}
            ls $PWD/*summary.txt > summaries.txt
        """

rule analyze_splicemutr:
    params:
        SUMMARY_TYPE=config["SUMMARY_TYPE"],
        NUM_SAMPLES=config["NUM_ALLELE_FILES"],
        SUMMARY_DIR=config["SUMMARY_DIR"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"]
    input:
        OUT_FILE_MHCNUGGETS=config["MHCNUGGETS_OUT"]+"/allele_files.txt",
        PROCESS_BINDAFF_FILES=config["PROCESS_BINDAFF_OUT"]+"/filenames.txt",
        EXTRACT_DATA_FILE=config["PROCESS_BINDAFF_OUT"]+"/summaries.txt",
        GENOTYPES=config["GENOTYPES_REFORMATTED"],
        SPLICE_DAT_FILE=config["SPLICE_DAT_FILE"],
        ANALYZE_SPLICEMUTR_OUT=config["ANALYZE_SPLICEMUTR_OUT"]
    output:
        ANALYZE_SPLICEMUTR_OUT_FILE=config["ANALYZE_SPLICEMUTR_OUT"]+"/filenames.txt"
    shell:
        """
            START=1
            for ((VAR=$START; VAR<={params.NUM_SAMPLES}; VAR++))
            do
                {params.SCRIPT_DIR}/analyze_splicemutr.py -g {input.GENOTYPES} -s {params.SUMMARY_DIR} -d {input.SPLICE_DAT_FILE} -o {input.ANALYZE_SPLICEMUTR_OUT} -t {params.SUMMARY_TYPE} -n $VAR
            done

            cd {input.ANALYZE_SPLICEMUTR_OUT}
            ls $PWD/*_splicemutr_kmers.txt > filenames.txt
        """

rule compile_kmer_counts:
    input:
        KMER_COUNTS_FILES=config["KMER_COUNTS_FILES"],
        SCRIPT_DIR=config["SPLICEMUTR_PYTHON"],
        KMER_COUNTS_OUT=config["KMER_COUNTS_OUT"]
    output:
        KMER_COUNTS_FILE=config["KMER_COUNTS_OUT"]+"/all_kmers_counts.txt"
    shell:
        """
        mkdir -p {input.KMER_COUNTS_OUT}

        {input.SCRIPT_DIR}/compile_kmer_counts.py -k {input.KMER_COUNTS_FILES} -o {input.KMER_COUNTS_OUT}
        """

rule create_junction_expression:
    input:
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        JUNC_FILES=config["JUNCFILES"],
        CREATE_JUNC_EXPRESSION_OUT=config["CREATE_JUNC_EXPRESSION_OUT"]
    output:
        CREATE_JUNC_EXPRESSION_FILE=config["CREATE_JUNC_EXPRESSION_OUT"]+"/junc_expr_combined_vst.rds"
    shell:
        """
        {input.SCRIPT_DIR}/create_junc_expr.R -f {input.JUNC_FILES} -o {input.CREATE_JUNC_EXPRESSION_OUT}
        """

rule calculate_gene_metric:
    input:
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        SPLICE_DAT_FILE=config["SPLICE_DAT_FILE"],
        KMER_COUNTS_FILE=config["KMER_COUNTS_FILE"],
        JUNC_EXPR_FILE=config["JUNC_EXPRESSSION_FILE"],
        CREATE_SPLICING_ANTIGENICITY_OUT=config["CREATE_SPLICING_ANTIGENICITY_OUT"]
    output:
        SPLICING_ANTIGENICITY_FILE=config["CREATE_SPLICING_ANTIGENICITY_OUT"]+"/filenames.txt"
    shell:
        """
        mkdir -p {input.CREATE_SPLICING_ANTIGENICITY_OUT}
        {input.SCRIPT_DIR}/calc_gene_metric_len_norm.R  -s {input.SPLICE_DAT_FILE} -k {input.KMER_COUNTS_FILE} -j {input.JUNC_EXPR_FILE} -o {input.CREATE_SPLICING_ANTIGENICITY_OUT}/CHOL
        cd {input.CREATE_SPLICING_ANTIGENICITY_OUT}
        ls $PWD/*_gene_metric_mean_len_norm_no_gene_norm_tumor.rds > filenames.txt
        """