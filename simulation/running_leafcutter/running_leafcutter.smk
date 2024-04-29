import os

configfile: "config.yaml"

if not os.path.exists(config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT"):
    os.makedirs(config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT")

rule all:
     input:
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/filenames.txt",
        RDATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/data.Rdata",
        JUNC_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt",
        GROUPS_FILE=config["GROUPS_FILE"],
        OUT_FILE_FINAL=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/splicemutr_introns.rds"

rule prepare_for_leafcutter:
    params:
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        SIMULATED_READS=os.getcwd()+"/"+config["SIMULATED_READS"]
    output:
        SJ_FILES=os.getcwd()+"/"+config["simulated_reads"]+"/bams/junc_files.txt"
    shell:
        """
        cd {params.SPLICEMUTR_SCRIPTS}
        chmod +x *
        cd {params.SIMULATED_READS}
        ls $PWD/*SJ.out.tab > junc_files.txt
        """

rule filtering_sj_files:
    input:
        SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/bams/junc_files.txt"
    output:
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/filenames.txt"
    params:
        FILTERED_SJ_FILES_OUT=directory(os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT"),
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        NUM_SAMPLES=config["NUM_SAMPLES"]
    shell:
        """
        START=1
        for ((VAR=$START; VAR<={params.NUM_SAMPLES}; VAR++))
        do
            SJ_FILE=$(sed -n ${{VAR}}p {input.SJ_FILES})
            {params.SCRIPT_DIR}/filter_juncs.R -o {params.FILTERED_SJ_FILES_OUT} -s $SJ_FILE
	    done

        cd {params.FILTERED_SJ_FILES_OUT}
        ls $PWD/*filt > filenames.txt
	    """

rule STAR_to_leaf:
    input:
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/filenames.txt"
    output:
        JUNC_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt"	
    params:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        FILTERED_SJ_FILES_OUT=directory(os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT"),
        FUNCTIONS=config["SPLICEMUTR_FUNCTIONS"],
        NUM_SAMPLES=config["NUM_SAMPLES"]
    shell:
        """
        START=1
        for ((VAR=$START; VAR<={params.NUM_SAMPLES}; VAR++))
        do
            SJ_FILE=$(sed -n ${{VAR}}p {input.FILTERED_SJ_FILES})
            {params.SCRIPT_DIR}/STAR_to_leaf.R -o {params.FILTERED_SJ_FILES_OUT} -s $SJ_FILE -f {params.FUNCTIONS}
        done
        
        cd {params.FILTERED_SJ_FILES_OUT}
        ls $PWD/*.junc > juncfiles.txt	
	    """ 

rule create_groups_file:
    output:
        GROUPS_FILE=os.getcwd()+"/"+config["GROUPS_FILE"]
    shell:
        """
        echo "sample_01.filt  1" >> {output.GROUPS_FILE}
        echo "sample_02.filt  1" >> {output.GROUPS_FILE}
        echo "sample_03.filt  1" >> {output.GROUPS_FILE}
        echo "sample_04.filt  1" >> {output.GROUPS_FILE}
        echo "sample_05.filt  2" >> {output.GROUPS_FILE}
        echo "sample_06.filt  2" >> {output.GROUPS_FILE}
        echo "sample_07.filt  2" >> {output.GROUPS_FILE}
        echo "sample_08.filt  2" >> {output.GROUPS_FILE}
        """

rule running_leafcutter:
    input:
        JUNC_DIR=directory(os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT"),
        JUNCFILE_FILENAMES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt",
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        LEAFCUTTER_SCRIPTS=os.getcwd()+"/"+config["LEAFCUTTER_SCRIPTS"],
        REF_DIR=config["REF_DIR"],
        LEAFVIZ_DIR=config["LEAFVIZ_DIR"],
        GROUPS_FILE=config["GROUPS_FILE"],
        LEAFCUTTER_PYTHON=config["LEAFCUTTER_PYTHON"],
    output:
        RDATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/data.Rdata"
    shell:
        """
        echo "leafcutter_cluster_regtools"
        python {input.LEAFCUTTER_PYTHON}/splicemutr_leafcutter_cluster_regtools.py -j {input.JUNCFILE_FILENAMES} -r {input.JUNC_DIR} -o data -l 500000

        echo "leafcutter_ds"
        {input.LEAFCUTTER_SCRIPTS}/leafcutter_ds.R -i 4 --num_threads 1 --exon_file={input.REF_DIR}/G039.exons.txt -o {input.JUNC_DIR}/leafcutter_ds {input.JUNC_DIR}/data_perind_numers.counts.gz {input.GROUPS_FILE}

        echo "prepare_results"
        {input.LEAFVIZ_DIR}/prepare_results.R -o {output.RDATA} -m {input.GROUPS_FILE} {input.JUNC_DIR}/data_perind_numers.counts.gz {input.JUNC_DIR}/leafcutter_ds_cluster_significance.txt {input.JUNC_DIR}/leafcutter_ds_effect_sizes.txt {input.REF_DIR}/G039
        """

rule save_introns:
    input:
        RDATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/data.Rdata",
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    output:
        OUT_FILE_FINAL=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["splicemutr"]+"/SJ_FILES_OUT/splicemutr_introns.rds"
    shell:
        """
        {input.SPLICEMUTR_SCRIPTS}/save_introns.R -i {input.RDATA} -o splicemutr
        """

