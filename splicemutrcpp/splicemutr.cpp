/*
The library function necessary for running the homolig input through cython from the python library. 
*/
#include "splicemutr.hpp"

const char *USAGE="usage: homolig [options] sequences.csv";
const char *HELP_INPUT="-h: help\n-i: seq_type\n-c: chains\n-m: metric\n-p: species";

enum ANN_TYPE check_annotation_type(const std::string &annotation_file){
    size_t annotation_type_gtf = annotation_file.find("gtf");
    size_t annotation_type_gff = annotation_file.find("gff3");
    if (std::string::npos!=annotation_type_gtf){
        return(GTF);
    } else if (std::string::npos != annotation_type_gff){
        return(GFF);
    } else {
        return(ERROR);
    }
};

std::string remove_leading_whitespace(std::string string_val){ // the whitespace after tab messes with casting otherwise
    if (string_val[0]==' '){
        string_val.substr(1,string_val.length());
    }
    return(string_val);
};

std::map<std::string,std::string> get_attributes(std::string attributes){
    std::map<std::string,std::string> line_info;
    std::stringstream line_stream(attributes);
    std::string line_split;
    while(std::getline(line_stream,line_split,';')){
        std::cout << line_split << std::endl;
        std::string line_split_filt = remove_leading_whitespace(line_split);
        std::stringstream subline_stream(line_split_filt);
        std::string subline_split_key,subline_split_value;
        std::getline(subline_stream,subline_split_key);
        std::getline(subline_stream,subline_split_value);
        line_info[subline_split_key]=subline_split_value;
    }
    return(line_info);
};

// std::string get_attribute(std::map<std::string,std::string> line_info,std::string attribute_to_get){
//     std::vector<std::string> line_info(GFF_NUM_ELEMS);
//     attributes = std::remove(attributes.begin(),attributes.end(),' ');
//     std::stringstream line_stream(attributes);
//     std::vector<std::string>::iterator line_elem=line_info.begin();
//     std::string line_split;
//     while(std::getline(line_stream,line_split,';')){
//         *line_elem = line_split;
//         line_elem++;
//     }
// }

GFF3 build_genes(std::ifstream &annotation_file_stream, enum ANN_TYPE ann, GFF3 GFF3_obj){ // this needs work for reading in reference for gene object
    if (ann == GFF){
        std::vector<std::string> line_info(GFF_NUM_ELEMS);
        std::string line;
        std::string line_split;
        bool IN_GENE_FLAG=false;
        while (std::getline(annotation_file_stream,line)){
            if (line.find("#") != std::string::npos){
                continue;
            }
            std::stringstream line_stream(line);
            std::vector<std::string>::iterator line_elem=line_info.begin();
            while(std::getline(line_stream,line_split,'\t')){
                std::string line_split_filt = remove_leading_whitespace(line_split);
                *line_elem = line_split_filt;
                line_elem++;
            }
            if (*(line_info.begin()+TYPE) == "gene" & IN_GENE_FLAG==false){
                IN_GENE_FLAG==true;
                int start_val = stoi(*(line_info.begin()+START));
                int end_val = stoi(*(line_info.begin()+END));
                std::string strand_val = *(line_info.begin()+STRAND);
                std::string chromosome = *(line_info.begin()+SEQID);
                std::vector<transcript> transcripts;
                std::string identity_val = *(line_info.begin()+ATTRIBUTES);
                gene gene_to_add = gene(start_val,end_val,strand_val,chromosome,transcripts,identity_val);
                GFF3_obj.add_gene(gene_to_add);
                std::cout << *(line_info.begin()+ATTRIBUTES) << std::endl;
            } else if (*(line_info.begin()+TYPE) == "gene" & IN_GENE_FLAG==true){
                std::cout << "inside" << std::endl;
                IN_GENE_FLAG=false;
                annotation_file_stream.seekg (-line_split.length(),std::ios_base::cur);
                continue;
            }
        }
        return(GFF3_obj);
    } else {
        return(GFF3_obj);
    }
};

bool splicemutr(const std::string &annotation_file,const std::string &splice_junction_file){
    enum ANN_TYPE annotation_type  = check_annotation_type(annotation_file);
    if (annotation_type==GTF){
        //do gtf priming stuff
    } else if (annotation_type==GFF){
        // do gff3 priming stuff
    } else {
        std::cerr << "annotation file does not have a .gtf or .gff extension" << std::endl;
        return(false);
    }
    std::ifstream annotation_file_stream(annotation_file);
    if (!annotation_file_stream.is_open()){
        std::cerr << "annotation file failed to open" << std::endl;
        return(false);
    }
    std::string line;
    std::string line_split;
    int line_num = 0;
    while (std::getline(annotation_file_stream,line)){ //scan past comments
        if (line.find("#") != std::string::npos){
            continue;
        } else {
            annotation_file_stream.seekg (-line.length(),std::ios_base::cur);
            break;
        }
    }
    //std::cout << line << std::endl;
    GFF3 GFF3_obj = GFF3();
    std::cout << "before build_genes()" << std::endl;
    GFF3_obj = build_genes(annotation_file_stream,annotation_type,GFF3_obj);
    std::cout << "after build_genes()" << std::endl;
    annotation_file_stream.close();
    return(true);
};