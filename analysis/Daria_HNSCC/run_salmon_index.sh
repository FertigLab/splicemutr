# running the salmon index method and saving it to head and neck data

TX=/media/theron/My_Passport/reference_genomes/SEQUENCES/GENCODE/gencode.v36.pc_transcripts.fa.gz

salmon index -t $TX --gencode -i gencodev35
