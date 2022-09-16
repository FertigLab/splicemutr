#ifndef SPLICEMUTR_HPP
#define SPLICEMUTR_HPP

#include "interval.hpp"
#include "genomic_annotation.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

extern const char *USAGE;
extern const char *HELP_INPUT;
const int GFF_NUM_ELEMS = 9;
enum ANN_TYPE { GTF, GFF, ERROR };
enum GFF_LINE_ELEMS { SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES };

ANN_TYPE check_annotation_type(const std::string &annotation_file);
std::string remove_first_character(std::string string_val);
bool build_gene(std::ifstream &annotation_file);
std::map<std::string,std::string> get_attributes(std::string attributes);
bool splicemutr(const std::string &annotation_file,const std::string &splice_junction_file);

#endif