/* 
A library for handling gff, gtf, and bed formatted files for splicemutr
Created: 20220830
Theron Palmer Jr.
*/

// GTF file format from: https://useast.ensembl.org/info/website/upload/GTF.html
// gtf file format from: https://useast.ensembl.org/info/website/upload/gff.html
// bed file format from: https://useast.ensembl.org/info/website/upload/bed.html

#include "genomic_annotation.hpp"

/********************************************************************************************************/
// The GTF class function definitions

GTF::GTF(void){};

GTF::GTF(std::string GTF_file_value){
    GTF_file = GTF_file_value;
};

bool GTF::check_genes(std::string gene_name){
    if (genes.count(gene_name)==0){
        return(true);
    } else {
        return(false);
    }
};

gene GTF::get_gene(std::string gene_name){
    try{
        gene gene_to_return = genes[gene_name];
        return(gene_to_return);
    } catch(std::out_of_range&){
        return(gene());
    }
};

bool GTF::add_gene(gene gene_to_add){
    std::string gene_to_add_identity = gene_to_add.get_identity();
    if (check_genes(gene_to_add_identity)){
        genes[gene_to_add_identity] = gene_to_add; // because is map
        return(true);
    } else {
        return(false);
    }
};

bool GTF::build_reference(void){ // this is where all of the text processing and use of the interval library happens
    return(true);
}; 

bool GTF::print(std::string output_file){ // this is so that we can print a modified GTF as a GTF file
    return(true);
}; 

size_t GTF::print_num_genes(void){
    return(genes.size());
};

std::vector<std::string> GTF::get_gene_names(void){
    std::vector<std::string> gene_names;
    for(std::map<std::string, gene>::iterator it = genes.begin(); it != genes.end(); ++it) {
        gene_names.push_back(it->first);
    }
    return(gene_names);
};

std::map<std::string,gene> GTF::get_genes(void){
    return(genes);
}

enum ANN_TYPE check_annotation_type(const std::string &annotation_file){
    size_t annotation_type_gtf = annotation_file.find("gtf");
    size_t annotation_type_gff = annotation_file.find("gff3");
    if (std::string::npos!=annotation_type_gtf){
        return(gtf);
    } else if (std::string::npos != annotation_type_gff){
        return(gff);
    } else {
        return(error);
    }
};

static inline std::string &remove_leading_whitespace(std::string &string_val) {
    string_val.erase(string_val.begin(), std::find_if(string_val.begin(), string_val.end(), [](int c) {return !std::isspace(c);}));
    return(string_val);
};

std::map<std::string,std::string> get_attributes(std::string attributes){ // this does not extract the key:value pair properly. It only extracts the key properly for the gene_id, or rather the first element. 
    std::map<std::string,std::string> line_info;
    std::stringstream line_stream(attributes);
    std::string line_split;
    while(std::getline(line_stream,line_split,';')){
        std::string line_split_filt = remove_leading_whitespace(line_split);
        std::stringstream subline_stream(line_split_filt);
        std::string subline_split_key,subline_split_value;
        std::string split_val;
        std::vector<std::string> split_vals;
        while(std::getline(subline_stream,split_val,' ')){
            split_vals.push_back(split_val);
        }
        subline_split_key = *(split_vals.begin());
        subline_split_value = *(split_vals.begin()+1);
        line_info[subline_split_key]=subline_split_value;
    }
    return(line_info);
};

std::vector<std::string> split_line(std::string& line_to_split){
    std::vector<std::string> line_info(GTF_NUM_ELEMS);
    std::string line_split;
    std::stringstream line_stream(line_to_split);
    std::vector<std::string>::iterator line_elem=line_info.begin();
    while(std::getline(line_stream,line_split,'\t')){
        std::string line_split_filt = remove_leading_whitespace(line_split);
        *line_elem = line_split_filt;
        line_elem++;
    }
    return(line_info);
};

void fill_gene_metadata(std::vector<std::string>& line_info, gene& gene_to_add){
    int start_val = stoi(*(line_info.begin()+START));
    int end_val = stoi(*(line_info.begin()+END));
    std::string strand_val = *(line_info.begin()+STRAND);
    std::string chromosome = *(line_info.begin()+SEQID);
    std::map<std::string,std::string> attributes_map = get_attributes(*(line_info.begin()+ATTRIBUTES));
    std::string identity_val = attributes_map["gene_id"];
    gene_to_add.set_start(start_val);
    gene_to_add.set_end(end_val);
    gene_to_add.set_strand(strand_val);
    gene_to_add.set_chromosome(chromosome);
    gene_to_add.set_identity(identity_val);
};

void fill_transcript_metadata(std::vector<std::string>& line_info, transcript& transcript_to_add){
    int start = stoi(*(line_info.begin()+START));
    int end = stoi(*(line_info.begin()+END));
    std::string strand = *(line_info.begin()+STRAND);
    std::string chromosome = *(line_info.begin()+SEQID);
    std::map<std::string,std::string> attributes_map = get_attributes(*(line_info.begin()+ATTRIBUTES));
    std::string identity_val = attributes_map["transcript_id"];
    transcript_to_add.set_start(start);
    transcript_to_add.set_end(end);
    transcript_to_add.set_strand(strand);
    transcript_to_add.set_chromosome(chromosome);
    transcript_to_add.set_identity(identity_val);
};

bool fill_transcript_element(std::vector<std::string>& line_info, transcript& transcript_to_fill){
    std::string element_type = *(line_info.begin()+TYPE);
    int start = stoi(*(line_info.begin()+START));
    int end = stoi(*(line_info.begin()+END));
    std::string strand = *(line_info.begin()+STRAND);
    std::string chromosome = *(line_info.begin()+SEQID);
    std::map<std::string,std::string> attributes_map = get_attributes(*(line_info.begin()+ATTRIBUTES));
    std::string identity_val = attributes_map["transcript_id"];

    if(element_type == "exon"){
        return(transcript_to_fill.add_exon(ginterval(start,end,chromosome,strand)));
    } else if (element_type == "CDS"){
        return(transcript_to_fill.add_cds(ginterval(start,end,chromosome,strand)));
    } else if (element_type == "five_prime_utr"){
        return(transcript_to_fill.add_utr5(ginterval(start,end,chromosome,strand)));
    } else if (element_type == "three_prime_utr"){
        return(transcript_to_fill.add_utr3(ginterval(start,end,chromosome,strand)));
    } else {
        //std::cerr << "failed to assign" << identity_val << ":" << element_type << ":" << chromosome << ":" << start << ":" << end << ":" << strand << std::endl;
        return(false);
    }
};

GTF build_gtf_ref(std::ifstream &annotation_file_stream, enum ANN_TYPE ann){ // this needs work for reading in reference for gene object
    GTF GTF_obj = GTF();
    if (ann == gtf){
        int num_genes = 0;
        std::string line;
        gene current_gene;
        transcript current_transcript;
        bool FIRST_GENE_FLAG=true;
        bool FIRST_TX_FLAG=true;
        while (std::getline(annotation_file_stream,line)){
            if (line.find("#") != std::string::npos){
                continue;
            } else {
                std::vector<std::string> line_info = split_line(line);
                if (*(line_info.begin()+TYPE) == "gene"){
                    if (FIRST_GENE_FLAG){
                        FIRST_GENE_FLAG=false;
                    } else {
                        if (!GTF_obj.add_gene(current_gene)){
                            std::cout << "failed to add " << current_gene.get_identity() << "with " << current_gene.get_transcript_identities().size() << " transcripts " << "to GTF_obj" << std::endl;
                        }
                    }
                    current_gene = gene();
                    fill_gene_metadata(line_info,current_gene);
                    if (num_genes==10){
                        break;
                    } else {
                        num_genes++;
                    }
                } else if (*(line_info.begin()+TYPE) == "transcript"){
                    if (FIRST_TX_FLAG){
                        FIRST_TX_FLAG=false;
                    } else {
                        std::cout << current_transcript.get_identity() << " : " << current_transcript.get_exons().get_gintervals().size() << " exons: " << current_transcript.get_cds().get_gintervals().size() << " cds: " << current_transcript.get_utr5().get_gintervals().size() << " 5'UTR: " << current_transcript.get_utr3().get_gintervals().size() << " 3'UTR" <<std::endl;
                        (current_gene).add_transcript(current_transcript);
                    }
                    current_transcript = transcript();
                    fill_transcript_metadata(line_info,current_transcript);
                } else {
                    fill_transcript_element(line_info,current_transcript);
                }
            }
        }
        return(GTF_obj);
    } else if (ann == gff){
        return(GTF_obj);
    } else {
        return(GTF_obj);
    }
};

void testing_GTF_ref(GTF GTF_obj){
    std::map<std::string,gene> genes = GTF_obj.get_genes();
    for(std::map<std::string, gene>::iterator it = genes.begin(); it != genes.end(); ++it) {
        std::string gene_name = it->first;
        gene current_gene = it->second;
        if (current_gene.is_empty()){
            std::cout << gene_name << ": is empty" << std::endl;
        } else {
            //std::cout << "gene->" << gene_name << ":" << "transcript->" << (current_gene.get_transcript_identities()).size() << std::endl;
        }
    }
    std::cout << GTF_obj.print_num_genes() << std::endl;
};