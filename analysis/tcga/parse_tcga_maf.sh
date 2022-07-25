# filtering mc3_maf for splicing factor genes

splice_factor_genes=/media/theron/My_Passport/TCGA_junctions/ext_dat/splicing_factor_genes.txt
maf_file=/media/theron/My_Passport/TCGA_junctions/mc3.v0.2.8.PUBLIC.maf
out_maf_file=/media/theron/My_Passport/TCGA_junctions/splice_factor.maf

sed -n 1p $maf_file > $out_maf_file
grep -F -f $splice_factor_genes $maf_file >> $out_maf_file
