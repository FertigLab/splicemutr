/* 
A library for handling gff4, gtf, and bed formatted files for splicemutr
Created: 20220830
Theron Palmer Jr.
*/

// gff3 file format from: https://useast.ensembl.org/info/website/upload/gff3.html
// gtf file format from: https://useast.ensembl.org/info/website/upload/gff.html
// bed file format from: https://useast.ensembl.org/info/website/upload/bed.html

#ifndef GENOMIC_ANNOTATION_HPP
#define GENOMIC_ANNOTATION_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "interval.hpp"

class GFF3 {
    protected:
        std::string GFF3_file;
        std::map<std::string, gene> genes;
    public:
        GFF3(void);
        GFF3(std::string GFF3_file_value);
        bool check_genes(std::string gene_name);
        gene get_gene(std::string gene_name);
        bool add_gene(gene gene_to_add);
        bool build_reference(void); // this is where all of the text processing and use of the interval library happens
        bool print(std::string output_file); // this is so that we can print a modified GFF3 as a GFF3 file
};

#endif