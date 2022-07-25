# running process_arcashla.py locally for analysis consideration

GENOTYPE_LOG_FILE=/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/DGay10-26144Aligned.out.bam_dir/DGay10-26144Aligned.genotype.log
GENOTYPE_JSON=/media/theron/My_Passport/head_and_neck_DARIA/data/splicemutr_05_26_2021/DGay10-26144Aligned.out.bam_dir/DGay10-26144Aligned.genotype.json

./process_arcashla.py -g $GENOTYPE_LOG_FILE -j $GENOTYPE_JSON