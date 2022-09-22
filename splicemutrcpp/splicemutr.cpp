/*
The library function necessary for running the homolig input through cython from the python library. 
*/
#include "splicemutr.hpp"

bool splicemutr(const std::string &annotation_file, const std::string &splice_junction_file){
    std::ifstream annotation_file_stream(annotation_file);
    if (!annotation_file_stream.is_open()){
        std::cerr << "annotation file failed to open" << std::endl;
        return(false);
    }
    enum ANN_TYPE annotation_type  = check_annotation_type(annotation_file);
    if (annotation_type==gtf){

        GTF GTF_obj = build_gtf_ref(annotation_file_stream,annotation_type);
        testing_GTF_ref(GTF_obj);

    } else if (annotation_type==gff){

    } else {
        std::cerr << "annotation file does not have a .gtf or .gff extension" << std::endl;
        return(false);
    }
    annotation_file_stream.close();
    return(true);
};