/* 
A library for handling gff4, gtf, and bed formatted files for splicemutr
Created: 20220830
Theron Palmer Jr.
*/

// GFF3 file format from: https://useast.ensembl.org/info/website/upload/GFF3.html
// gtf file format from: https://useast.ensembl.org/info/website/upload/gff.html
// bed file format from: https://useast.ensembl.org/info/website/upload/bed.html

#ifndef GENOMIC_ANNOTATION_HPP
#define GENOMIC_ANNOTATION_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "interval.hpp"
#include <sstream>
#include <algorithm>

class GTF {
    protected:
        std::string GTF_file;
        std::map<std::string, gene> genes;
    public:
        GTF(void);
        GTF(std::string GTF_file_value);
        bool check_genes(std::string gene_name);
        gene get_gene(std::string gene_name);
        bool add_gene(gene gene_to_add);
        bool build_reference(void); // this is where all of the text processing and use of the interval library happens
        bool print(std::string output_file); // this is so that we can print a modified GTF as a GTF file
        size_t print_num_genes(void);
        std::vector<std::string> get_gene_names(void);
        std::map<std::string,gene> get_genes(void);
};

const int GTF_NUM_ELEMS = 9;
enum ANN_TYPE {gtf,gff,error};
enum GFF_LINE_ELEMS { SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES };

enum ANN_TYPE check_annotation_type(const std::string &annotation_file);
static inline std::string &remove_leading_whitespace(std::string &string_val);
std::map<std::string,std::string> get_attributes(std::string attributes);
std::vector<std::string> split_line(std::string& line_to_split);
void fill_gene_metadata(std::vector<std::string>& line_info, gene& gene_to_add);
void fill_transcript_metadata(std::vector<std::string>& line_info, transcript& transcript_to_add);
bool fill_transcript_element(std::vector<std::string>& line_info, transcript& transcript_to_fill);
GTF build_gtf_ref(std::ifstream &annotation_file_stream, enum ANN_TYPE ann);
void testing_GTF_ref(GTF GTF_obj);

#endif