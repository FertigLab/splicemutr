# SpliceMutr: Calculating splicing-derived neoantigen burden from splice-junctions affecting protein-coding transcripts in RNA-seq data.

------------------------------------------------------------------------

SpliceMutr is a tool used for the analysis of splicing-derived antigen burden. SpliceMutr uses differential intron usage analysis through LeafCutter or LeafCutterMD to first determine differential usage of splice-junctions between sets of RNA-seq samples. LeafCutter uses a Dirichlet-Multinomial generalized linear model to fit splice-junction counts within clusters of splice-junctions grouped due to overlapping splice site usage between pairs of sample types (tumor and normal or pre and post treatment in the case of the original SpliceMutr paper). Two separate Dirichlet-Multinomial models are fit to the splice-junction counts per splice-junction cluster based on the sample type stratification and a likelihood ratio test is used to evaluate differential usage within the splice-junction cluster as a whole. Each individual differentially-used splice junction per cluster is then used to modify the reference transcriptome. The modified transcripts are searched for canonical and/or modified open reading frames, translated based on the best identified open reading frame, and kmerized into 9-mers. During kmerization, a map between the peptide kmer and modified protein of interest is generated for searching after MHC-peptide binding prediction. MHC-peptide binding prediction is carried out using the arcasHLA generated HLA genotype per sample analyzed, the set of generated peptide 9-mers, and MHCnuggets or NetMHCpan. Peptide 9-mers found to bind to any of the HLA molecules per sample are then mapped back to their transcript of origin using the kmer-transcript map and peptide binding summaries are generated that contain the number of unique binding kmers per transcript and per sample. If differential intron usage analysis is performed, the number of unique and sample-type specific kmers are generated by filtering out all kmers shared between sample types within splice-junction clusters. The number of these unique peptide kmers per splice-junction-modified transcript is used to generate a splicing antigenicity score per sample and per gene that is higher if the splice-junctions associated with a gene produce more peptide kmers able to bind to an MHC complex and have higher RNA-seq expression. This score can then be associated to other quantitative measures.

------------------------------------------------------------------------

# Installation

Requirements

[STAR](https://github.com/alexdobin/STAR)\
[arcasHLA](https://anaconda.org/bioconda/arcas-hla)\
[MHCnuggets](https://karchinlab.org/apps/appMHCnuggets.html) or [NetMHCpan](https://github.com/tzina97/netMHCpan/tree/main)\
[LeafCutter](https://davidaknowles.github.io/leafcutter/)\
R \>= 4.0.2\
python \>= 3.6.10\

To complete installation install the splicemutr package within your installation of R. This will install the packages necessary and prompt you to install the associated bioconductor packages.

``` r
install_github("FertigLab/splicemute")
```

To build the SpliceMutr environment used in the SpliceMutr paper, install [miniconda](https://docs.conda.io/projects/miniconda/en/latest/), create a conda environment, and build the conda package contained within splicemutr.yml

``` bash
conda env create -f environment.yml
```

------------------------------------------------------------------------

SpliceMutr is licensed under the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html)

If you find SpliceMutr useful, please cite as follows:

Palmer, Theron, et al. "SpliceMutr Enables Pan-Cancer Analysis of Splicing-Derived Neoantigen Burden in Tumors." bioRxiv, 2023, <https://doi.org/10.1101/2023.05.26.542165>.

Examples of how the SpliceMutr pipeline scripts have been used for the analysis associated with the manuscript "SpliceMutr enables pan-cancer analysis of splicing-derived neoantigen burden in tumors" can be found in the folder "analysis" under the TCGA and CM038 headings and within the pipeline section. There is a jupyter notebook explaining usage of the appropriate scripts.

------------------------------------------------------------------------

# Usage

\*\*\* the snakemake file contained within this directory is for reference. It will not run a sample without modification and this pipeline will not run to completion without your input primarily through formatting of intermediate files and necessary files for varying pipeline components.

1)  This step takes in a GTF file and creates the sqlite database necessary for transcript formation in form_transcripts.R. This is covered under rule make_txdb in the snakemake file within this directory.

``` bash
make_txdb.R
-o output_file.sqlite
-g GTF
```

output_file.sqlite - the sqlite file to write\
GTF - the GTF file to convert into the txdb database (sqlite)\

2)  Prepare the BSgenome references used for transcript formation in form_transcripts.R. This is covered first by rule prep_bsgenome_reference in the snakemake file within this directory which covers converting the BSgenome fasta file to a 2bit file.Building the BSgenome package is done using the following directions: <https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf>

3)  Prepare the LeafCutter references. LeafCutter reference preparation is covered in the rule prepare_leafcutter_references in the snakemake file within this directory. The directions outlined at [LeafCutter](https://davidaknowles.github.io/leafcutter/) are what are followed to create the LeafCutter references.

4)  Align your bulk RNA-seq reads to the appropriate genome of choice using STAR. Use the directions outlined at [LeafCutter](https://davidaknowles.github.io/leafcutter/) for building your STAR index as well as alignment of your reads to the STAR index you built. For the SpliceMutr paper analysis, we used prebuilt indeces that are no longer accessible at the associated LeafCutter link in the directions.

5)  This step takes the STAR sj.out.tab files output from STAR alignment and formats them into the .junc file format that LeafCutter accepts. The rule convert_STAR_sj.out.tab_to_leafcutter\_.junc in the snakemake within this directory covers this protocol.

``` bash
STAR_to_leaf.R
-o output_directory
-s star_file
-f splicemutr_functions
```

output_directory - the directory to write the .junc file to\
star_file - the sj.out.tab file to convert into the .junc file\
splicemutr_functions - the splicemutr functions.R file to load into the script\

6)  Run LeafCutter on the .junc files generated for your experiment. You will need to format the groups and junctions files for your experiments. The groups file is a tab-delimited file with junc file name in one column and binary experimental grouping in the second column. The first experimental group in the groups file must be your reference group for the analysis. Further details about running LeafCutter are explained at [LeafCutter](https://davidaknowles.github.io/leafcutter/) and are covered in the rule running_leafcutter in the snakemake within this directory.

7)  This step takes the LeafCutter LeafViz output and saves the introns in a format suitable for input into SpliceMutr format_transcripts.R. This is covered in rule save_introns in the associated snakemake file within this directory.

``` bash
save_introns.R
-i intron_file.rds
-o out_prefix
```

intron_file.rds - the intron file output from LeafCutter LeafViz differential intron usage analysis out_prefix - the outprefix for the introns .rds files output from save_introns.R

7.1) I often use split_introns.R for analysis. This script takes similar inputs as save_introns.R, but it also takes a split_num argument for splitting the introns data.rds file output from LeafCutter LeafViz into split_num chunked intron.rds files.

``` bash
split_introns.R
-i intron_file.rds
-o out_prefix
-s split_num
```

intron_file.rds - the intron file output from LeafCutter LeafViz differential intron usage analysis out_prefix - the out prefix for the introns .rds files output from split_introns.R split_num - the number of introns to chunk the LeafCutter Leafvix .rds file into

8)  This step modifies transcripts using the identified differentially-used introns from LeafCutter. The rule form_transcripts in the snakemake within this directory covers usage.

``` bash
form_transcripts.R
-o out_prefix
-t txdb_file
-j junc_file
-b bsgenome_name
-m chr_map
-f funcs
```

out_prefix - the output prefix for the\
txdb_file - the sqlite database created in the rule make_txdb and created from make_txdb.R\
junc_file - the introns.rds file created from the rule save_introns (save_introns.R) or one of the files from split_introns.R if batching your job\
bsgenome_name - the name of the bsgenome that you created from the reference fasta file.\
chr_map - logical indicating whether you want to replace "chr" in the junc_file (introns.rds) with "". This is sometimes necessary if your gtf or fasta references do not match regarding this chromosome naming convention.\
funcs - the path to the SpliceMutr functions.R file\

9)  This step calcualtes the LGC coding potential of the transcripts generated from form_transcripts. The rule calculate_coding_potential in the snakemake within this directory covers usage.

``` bash
calc_coding_potential.R
-s splicemutr_data
-t transcript_fasta
-o out_dir
-f funcs
```

splicemutr_data - the \*\_splicemutr.rds file containing the splicemutr metadata output from form_transcripts.R\
transcript_fasta - the fasta output associated with the splicemutr metadata ouput from form_transcripts.R\
out_dir - the output directory to write the coding-potential-corrected SpliceMutr metadata ouptut from form_transcripts.R\
funcs - the path to the SpliceMutr functions.R file\

10) This step combines all batched (if batched) SpliceMutr metadata files and generates a file containing the set of proteins per modified transcript and a file containing the splicemutr metadata subset by those modified transripts that can form proteins. The rule combine_splicemutr within the snakemake within this directory covers usage.

``` bash
combine_splicemutr.R
-o out_dir
-s splice_files
```

out_dir - The output directory that the combined SpliceMutr files are written to\
splice_files - The splicemutr files generated by\
form_transcripts.R in a .txt file where each line is a file, including its directory\

11) This step processes the proteins.txt file from the previous step and generates a kmer-transcript location (row within the peptides.txt file) map. It also generates a peps_kmer_length.txt file that contains the unique set of k-mers from the kmer-transcript map. The rule process_peptides gives further details about usage.

``` bash
process_peptides.py
-p peptides.txt
-o output_directory
-k kmer_length
```

peptides.txt - the peptides.txt file output from combine_splicemutr.R\
output_directory - the directory to write the kmer-transcript map dictionary and the peps_9.txt file\
kmer_length - the kmer length to use for the kmerization of the peptides in petides.txt

12) This step runs arcasHLA on the STAR generated .bam files. arcasHLA performs HLA genotyping on bulk RNA-seq data. The output of arcasHLA will need to be formatted using the genotype.txt format outlined [here](). The rule, run_arcasHLA in the snakemake within this directory covers running arcasHLA.

13) This step runs MHcnuggets on the HLA genotypes, or superalleles associated with the HLA genotypes, output from arcasHLA and formatted as a file with one allele per line using the genotypes.txt HLA format and the peps_kmer_length.txt file generated from combine_splicemutr.R.

``` bash
runMHCnuggets.py
-t type
-k kmers
-m mhc
-o output
```

type - the HLA allele type "I" or "II" kmers - the peps_kmer_length.txt file containing the set of unique peptide kmers mhc - the mhc allele file with one allele per line output - the output directory to write the mhcnuggets results to

14) This step processes the binding affinity prediction results output from MHCnuggets and filters the kmers for those that have predicted IC50 scores \<=500 nM. The rule process_bindaffinity details usage for this step.

``` bash
process_bindaff.py
-b binders
-p pickle_dir
-o out
-k kmer_length
```

binders - a MHCnuggets .txt output file pickle_dir - the directory that the pickle file is written to out - the output directory to write the filtered binders file to kmer_length - the kmer_length to process

15) This step extracts the binding location data from the pickle and from the associated allele kmers and generates a summary file per allele that contains location and binding kmer information for the specific allele. The rule extract_data in the snakemake within this directory details usage for this step.

``` bash
extract_data.py
-a hla_allele
-p pickle_directory
-b kmer_beginning_range
-e kmer_end_range
```

hla_allele - the hla allele in genotypes.txt hla format to extract binding kmers for\
pickle_directory - the location of the kmer-transcript location map file\
kmer_beginning_range - the beginning kmer size to analyse when generating the HLA summary file. Typically set to 9\
kmer_end_range - the end kmer size to analyse when generating the HLA summary file. Only the number directly before this kmer size is the end size analyzed. Typically set to 10.

16) This step generates a kmer counts file per sample analyzed that is later compiled into a counts per sample file used during splicing antigenicity calculations. The rule analyze_splicemutr within the snakemake within this directory details usage for this step.

``` bash
-g genotypes_file
-s summary_dir
-d splice_dat
-o output_dir
-t summary_type
-n sample_num
```

genotypes_file - the HLA genotypes per sample in the genotypes.txt file format\
summary_dir - the directory containing the summary files for each HLA in the genotypes.txt file\
splice_dat - the data_splicemutr_all_pep.txt file output from combine_splicemutr.R\
output_dir - the output directory to write the associated sample kmer counts file to\
summary_type - 'IC50' or 'perc'. IC50 reads in a summary file with the format \*\_tx_dict_summary.txt while 'perc' reads in a summary file with the format \*\_tx_dict_summary_perc.txt.
sample_num - an integer indicating the row number of/sample to process in the genotype.txt file\

17) This step takes the kmer files generated by analyze_splicemutr.py run batched for all samples analyzed and compiles them into a single file where each column is a sample and each row corresponds to a row in the splicemutr *all_pep.rds metadata file output from combine_splicemutr.R. The rule compile_kmer_counts in the snakemake within this directory details usage for this step. 

``` bash
compile_kmer_counts.py
-k kmers_files
-o output_dir
```

kmers_files - a file containing the kmer counts files per sample, one each line with their directories\
output_dir - the output directory to write the all_kmers_counts.txt file to.\

18) This step uses the .junc files generated using STAR_to_leaf.R and creates a variance-stabilized junction expression matrix with junctions as rows and samples as columns. The rule create_junction_expression in the snakemake within this directory details usage for this step. 

``` bash
create_junc_expression.R
-j junc_dir
-f junc_files
-o output_directory
```

junc_dir - the junctions file directory\
junc_files - the junction files and their directories, one per line\
output_directory - the output directory to write the junction expression to\


19) This step calculates the splicing antigenicity per gene and sample. This script also generates the average coding potential per gene and sample. The rule calculate_gene_metric in the snakemake within this directory details usage for this step. 

``` bash
calculate_gene_metric_len_norm.R
-s splice_dat_file
-k kmer_counts
-j junc_expr_file
-o output_directory
```

splice_dat_file - the *all_pep.rds SpliceMutr file output from combine_splicemutr.R\
kmer_counts - the all_kmers_counts.txt file output from compile_kmer_counts.py\
junc_expr_file - the file containing the combined junction expression file output from create_junc_expression\
output_directory - the output directory to write the SpliceMutr splicing antigenicity calculation files to\
