import os

configfile: "config.yaml"

if not os.path.exists(config["FORMED_TRANSCRIPTS_DIR"]):
    os.mkdir(config["FORMED_TRANSCRIPTS_DIR"])
if not os.path.exists(config["COMBINE_SPLICEMUTR_OUT"]):
    os.mkdir(config["COMBINE_SPLICEMUTR_OUT"])
if not os.path.exists(config["PROCESS_PEPTIDES_OUT"]):
    os.mkdir(config["PROCESS_PEPTIDES_OUT"])
if not os.path.exists(config["GENOTYPES_DIR"]):
    os.mkdir(config["GENOTYPES_DIR"])

#rule all:
#    input:
        #FORMED_TRANSCRIPTS=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr.rds"
        #FORMED_TRANSCRIPTS_CP=config["FORMED_TRANSCRIPTS_DIR"]+"/CHOL_introns_data_splicemutr_cp_corrected.rds"
        #OUTPUT_FILE=config["COMBINE_SPLICEMUTR_OUT"]+"/data_splicemutr_all_pep.txt"
        #OUT_FILE=config["PROCESS_PEPTIDES_OUT"]+"/peps_9.txt"

'''   
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
        SPLICE_FILE=config["SPLICEMUTR_FILE"],
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
'''

rule run_arcasHLA:
  input:
    GENOTYPES_DIR=config["GENOTYPES_DIR"],
    FILENAMES_FILE=config["BAMFILES"],
  output:
    
  shell:
    """
    conda activate miniconda3/envs/arcashla
    START=1
    END=2
    for VAR in {{$START..$END}}
    do
      FILE=$(sed -n ${{VAR}}p {input.FILENAMES_FILE})
      FILE_BASE=$(basename {input.FILE})
      FILE_DIR={input.GENOTYPES_DIR}/${{FILE_BASE}}_dir
      mkdir $FILE_DIR

      # sort bam file
      samtools sort -o ${{FILE}}.sorted $FILE

      arcasHLA extract ${{FILE}} -o $FILE_DIR -v

      cd $FILE_DIR

      FASTQ1=$(ls *.extracted.1*)
      FASTQ2=$(ls *.extracted.2*)
      arcasHLA genotype $FASTQ1 $FASTQ2 -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRA,DRB1 -o $FILE_DIR -v
    done
    """