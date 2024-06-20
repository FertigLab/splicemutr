import os

configfile: "config.yaml"

if not os.path.exists(config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT"):
    os.makedirs(config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT")

rule all:
     input:
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/filenames.txt",
        DATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/data_pVals.txt",
        JUNC_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt",
        OUT_FILE_FINAL=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/comparison_juncs_linear.rds"

rule prepare_for_leafcutter:
    params:
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        SIMULATED_READS=os.getcwd()+"/"+config["SIMULATED_READS"]+"/bams"
    output:
        SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/bams/junc_files.txt"
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
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/filenames.txt"
    params:
        FILTERED_SJ_FILES_OUT=directory(os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT"),
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
        FILTERED_SJ_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/filenames.txt"
    output:
        JUNC_FILES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt"	
    params:
        SCRIPT_DIR=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"],
        FILTERED_SJ_FILES_OUT=directory(os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT"),
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

rule running_leafcutter:
    input:
        JUNC_DIR=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT",
        JUNCFILE_FILENAMES=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/juncfiles.txt",
        SPLICEMUTR_PYTHON=os.getcwd()+"/"+config["SPLICEMUTR_PYTHON"],
        LEAFCUTTER_SCRIPTS=os.getcwd()+"/"+config["LEAFCUTTER_SCRIPTS"],
        LEAFCUTTER_PYTHON=os.getcwd()+"/"+config["LEAFCUTTER_PYTHON"]
    output:
        DATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/data_pVals.txt"
    shell:
        """
        echo "leafcutter_cluster_regtools"
        python {input.SPLICEMUTR_PYTHON}/splicemutr_leafcutter_cluster_regtools.py -j {input.JUNCFILE_FILENAMES} -r {input.JUNC_DIR} -o data -l 500000

        cd {input.JUNC_DIR}

        echo "leafcutter_outlier"
        {input.LEAFCUTTER_SCRIPTS}/leafcutterMD.R --num_threads 1 -o data data_perind_numers.counts.gz

        """

rule save_introns:
    input:
        DATA=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/data_pVals.txt",
        SPLICEMUTR_SCRIPTS=os.getcwd()+"/"+config["SPLICEMUTR_SCRIPTS"]
    output:
        OUT_FILE_FINAL=os.getcwd()+"/"+config["SIMULATED_READS"]+"/"+config["SPLICEMUTR"]+"/SJ_FILES_OUT/comparison_juncs_linear.rds"
    shell:
        """
        chmod +x {input.SPLICEMUTR_SCRIPTS}/*
        {input.SPLICEMUTR_SCRIPTS}/create_outlier_linear_juncs.R -o {input.DATA}
        """
