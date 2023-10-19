configfile: "config.yaml"

rule convert_STAR_sj.out.tab_to_leafcutter_.junc:
  input:
    SJ_FILES=BAMS/SJ_files.txt
    SCRIPT_DIR=splicemute/scripts
    FUNCTIONS=/splicemute/R/functions.R
    NUM_SJ_FILES=8
  output:
    OUT_DIR="junc_files"
  shell:
    "conda activate miniconda3/envs/splicemutr

    START=1
    for VAR in {$START..{input.NUM_SJ_FILES}}
    do
      STAR_JUNCFILE=$(sed -n ${VAR}p {input.SJ_FILES})
      {input.SCRIPT_DIR}/STAR_to_leaf.R -f {input.FUNCTIONS} -o {input.OUT_DIR} -s {input.STAR_JUNCFILE}
    done

    cd {output.OUT_DIR}
    ls $PWD/*.junc > filenames.txt

    conda deactivate"

rule running_leafcutter:
  input:
      JUNC_DIR=/junc_files
      LEAF_SCRIPTS=/leafcutter/scripts
      REF_DIR=references/leafcutter_annotations
      LEAFVIZ_DIR=/leafcutter/leafviz
      GROUPS_FILE=/groups_file.txt
  shell:
    conda activate miniconda3/envs/splicemutr

    echo "leafcutter_cluster_regtools"
    python2 {input.LEAF_SCRIPTS}/splicemutr_leafcutter_cluster_regtools.py -j {input.JUNC_DIR}/junc_file.txt -r {input.JUNC_DIR} -o data -l 500000

    echo "leafcutter_ds"
    {input.LEAF_SCRIPTS}/leafcutter_ds.R --num_threads 1 --exon_file={input.REF_DIR}/G026.exons.txt.gz -o {input.JUNC_DIR}/leafcutter_ds {input.JUNC_DIR}/data_perind_numers.counts.gz {input.GROUPS_FILE}

    echo "prepare_results"
    {input.LEAFVIZ_DIR}/prepare_results.R -o {input.JUNC_DIR}/data.Rdata -m {input.GROUPS_FILE} {input.JUNC_DIR}/data_perind_numers.counts.gz {input.JUNC_DIR}/leafcutter_ds_cluster_significance.txt {input.JUNC_DIR}/leafcutter_ds_effect_sizes.txt {input.REF_DIR}/G026

    conda deactivate"

rule save_introns:
  input:
    INTRON_FILE=/junc_files/data.Rdata
    SPLIT_NUM=5000
    SCRIPT_DIR=splicemute/scripts
  output:
    OUT_DIR=/introns
  shell:
    "conda activate miniconda3/envs/splicemutr

    {input.SCRIPT_DIR}=/users/tpalmer/splicemute/scripts

    {input.SCRIPT_DIR}/save_introns.R -i {input.INTRON_FILE} -o {output.OUT_DIR}

   conda deactivate"