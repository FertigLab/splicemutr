
Each .sh file listed below outlines the usage of the associated .R or .py script in the "scripts" or "inst" subdirectories specifically for producing the data used in the manuscript. \

# building the necessary reference

See the "reference" subdirectory\

# Running SpliceMutr for the TCGA data

runLeafcutter.sh - Running LeafCutter differential splicing analysis\
recount3_tcga_juncs_init.sh - Downloading the recount3 splice junction data\
TCGA_gene_expr_per_sample.sh - Downloading the recount3 gene expression per TCGA sample\
create_junc_expr.sh - Creating the variance-stabilized and non-variance stabilized splice junction expression needed for input into SpliceMutr splicing antigenicity calculations \
TCGA_create_genotypes.sh - Creating the genotype information used during MHC:peptide binding affinity calculations. The genotype information come from https://gdc.cancer.gov/about-data/publications/panimmune and are found in the GDC-controlled file "OptiTypeCallsHLA_20171207.tsv" \
form_transcripts.sh - Performing the SpliceMutr transcript formation \
calc_coding_potential.sh - Recalculating the coding potential. This was necessary for the manuscript because I originally programmed the calculation incorrectly. \"form_transcripts.R" calculates the coding potential correctly now. \
combine_splicemutr.sh - Typically, SpliceMutr transcript formation is batched. This combines the output of each of those batched jobs. \
process_peptides.sh - Performs kmerization and tracks the protein of origin for the proteins output from SpliceMutr transcript formation. This ensures there is no need to realign kmers to the formed proteins after MHC:peptide binding affinity prediction. The .pickle files created help track the protein of origin for each kmer.  \
runMHCnuggets.sh - Performs MHCnuggets MHC:peptide binding affinity predictions on the kmerized peptides using the collapsed (all unique alleles for all samples analyzed) genotype.\
process_bindaff.sh - Filters the output from MHCnuggets MHC:Peptide binding affinity predictions for binding kmers. \
extract_data_large.sh - Extracts the set of binding kmers per allele and maps binders to samples based on sample genotype.\
analyze_splicemutr/prep_*.sh - Prepares the directory that will be written to during analyze_splicemutr_*.sh\
analyze_splicemutr/analyze_splicemutr_*.sh - Creates the set of files necessary for calculating the splicing antigenicity per TCGA tumor subtype.\
calc_gene_metric.sh - Calculates the splicing antigenicity per TCGA tumor subtype, sample, and gene\
splicemutr/analysis/TCGA/analysis_Theron_10132022.Rmd - Performs the analysis and generates the figures found in the manuscript. \
