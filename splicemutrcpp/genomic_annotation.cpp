/* 
A library for handling gff, gtf, and bed formatted files for splicemutr
Created: 20220830
Theron Palmer Jr.
*/

// gff3 file format from: https://useast.ensembl.org/info/website/upload/gff3.html
// gtf file format from: https://useast.ensembl.org/info/website/upload/gff.html
// bed file format from: https://useast.ensembl.org/info/website/upload/bed.html

#include "genomic_annotation.hpp"

/********************************************************************************************************/
// The gene class function definitions

GFF3::GFF3(void){};

GFF3::GFF3(std::string GFF3_file_value){
    GFF3_file = GFF3_file_value;
};

bool GFF3::check_genes(std::string gene_name){
    try {
        genes[gene_name];
    }   catch(std::out_of_range&){
        return(false);
    }
    return(true);
};

gene GFF3::get_gene(std::string gene_name){
    try{
        gene gene_to_return = genes[gene_name];
        return(gene_to_return);
    } catch(std::out_of_range&){
        return(gene());
    }
};

bool GFF3::add_gene(gene gene_to_add){
    std::string gene_to_add_identity = gene_to_add.get_identity();
    if (!check_genes(gene_to_add_identity)){
        genes[gene_to_add_identity] = gene_to_add; // because is map
        return(true);
    } else {
        return(false);
    }
};

bool GFF3::build_reference(void){ // this is where all of the text processing and use of the interval library happens
    return(true);
}; 

bool GFF3::print(std::string output_file){ // this is so that we can print a modified GFF3 as a GFF3 file
    return(true);
}; 