import os

configfile: "config.yaml"

if not os.path.exists(config["FILTERED_SJ_FILES_OUT"]):
    os.mkdir(config["FILTERED_SJ_FILES_OUT"])

rule all:
     input:
        RDATA=config["FILTERED_SJ_FILES_OUT"]+"/data.Rdata",
        SJ_FILES=config["FILT_SJ_FILES"],
	JUNC_FILES=config["JUNCFILE_FILENAMES"],
        OUT_FILE_FINAL=config["FILTERED_SJ_FILES_OUT"]+"/splicemutr_introns.rds"
	
rule filtering_sj_files:
    input:
        SJ_FILES=config["SJ_FILES"]
    output:
        FILTERED_SJ_FILES=config["FILT_SJ_FILES"]
    params:
        FILTERED_SJ_FILES_OUT=config["FILTERED_SJ_FILES_OUT"],
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
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
        FILTERED_SJ_FILES=config["FILT_SJ_FILES"]
    output:
        JUNC_FILES=config["JUNCFILE_FILENAMES"]	
    params:
        SCRIPT_DIR=config["SPLICEMUTR_SCRIPTS"],
        FILTERED_SJ_FILES_OUT=config["FILTERED_SJ_FILES_OUT"],
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
        JUNC_DIR=config["FILTERED_SJ_FILES_OUT"],
        JUNCFILE_FILENAMES=config["JUNCFILE_FILENAMES"],
        LEAFCUTTER_SCRIPTS=config["LEAFCUTTER_SCRIPTS"],
        REF_DIR=config["REF_DIR"],
        LEAFVIZ_DIR=config["LEAFVIZ_DIR"],
        GROUPS_FILE=config["GROUPS_FILE"],
        LEAFCUTTER_PYTHON=config["LEAFCUTTER_PYTHON"],
    output:
        RDATA=config["FILTERED_SJ_FILES_OUT"]+"/data.Rdata"
    shell:
        """
        echo "leafcutter_cluster_regtools"
        python {input.LEAFCUTTER_PYTHON}/splicemutr_leafcutter_cluster_regtools.py -j {input.JUNCFILE_FILENAMES} -r {input.JUNC_DIR} -o data -l 500000

        echo "leafcutter_outlier"
        $LEAF_SCRIPTS/leafcutterMD.R --num_threads 1 -o sample_01.filt sample_01.filt_perind_numers.counts.gz

        """