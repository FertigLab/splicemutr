Each .sh file listed below outlines the usage of the associated .R or .py script in the "scripts" or "inst" subdirectories specifically for producing the data used in the manuscript.

# Building the necessary reference

See the "reference" subdirectory

# Running SpliceMutr for the CM038 data

run_arcashla.sh - Running arcasHLA on the CM038 RNA-seq data to generate the per-sample genotype.

fastq2bams_valsamo.sh - Aligning CM038 RNA-seq data to the reference genome.

filter_juncs.sh - Filtering out non-canonical splice sites from the SJ.out.tab files.

star_to_leaf.sh - Converting the filtered STAR Sj.out.tab files to .junc files for input into LeafCutterMD.

splicemutr_leafcutter_valsamo.sh - Running LeafCutterMD for the CM038 trial RNA-seq data.

form_transcripts.sh - Performing SpliceMutr transcript formation.

calc_coding_potential.sh - Recalulcating the coding potential. I failed in programming the correct coding potential calculation during my first attempt. form_transcripts.R now calculates the coding potential properly. 

combine_splicemutr_cp.sh - Combining the batched output of SpliceMutr form_transcripts.R.

process_peptides.sh - Processing the proteins output from form_transcripts.R and generating a kmer map so that binding kmers do not need to be aligned in order to find their protein origin. 

runMHCnuggets.sh - Running MHCnuggets on the unique set of KMERS generated from process_peptides.sh per associated class 1 MHC allele for the CM038 samples.

process_bindaff.sh - Filtering out non-binders per HLA allele.

extract_data.sh - Creating per-allele binder counts for each SpliceMutr formed transcript.

run_featurecounts.sh - Running featurecounts on CM038 alignment data.

combine_featurecounts.sh - Combining featurecounts data run on the CM038 data.

create_junc_expression.sh - Combining .junc expression per sample into one object for use with the splicing antigenicity calcualtions.

valsamo_analyze_splicemutr.sh - Calculating the binding kmer counts per SpliceMutr-formed transcript and per sample. 

valsamo_create_comparisons_cp.sh - Using CM038 and SpliceMutr metadata to create the files necessary for the splicing antigenicity calculation per comparison.

updated_calc_gene_metric_len_norm.sh - Calculating the splicing antigenicity using the output from valsamo_create_comparisons_cp.sh.
 
