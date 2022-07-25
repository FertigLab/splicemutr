# The purpose of this script is to run process_bindaff.R on the splicemutr formatted netmhc data

SCRIPTS=/media/theron/My_Passport/splicemutr
KMER_DIR=/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/kmers
SCORE_DIR=/media/theron/My_Passport/head_and_neck_DARIA/data/leafcutter_10_21_2020/leafcutter_introns_split/predictions_1
KMER_LEN=9
LINE_NUM=5

$SCRIPTS/process_bindaff.R -o $KMER_DIR -i $SCORE_DIR -k $KMER_LEN -m $LINE_NUM
