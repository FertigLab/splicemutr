/* 
a library for handling genomic intervals and sets of genomic intervals (ie transcripts and genes)
Created: 20220824
Theron Palmer Jr.
*/

#include "interval.hpp"

/********************************************************************************************************/
// The interval class function definitions

interval::interval(){}; // constructor

interval::interval(int start_value,int end_value){ // constructor
    start = start_value;
    end = end_value;
};

void interval::set_start(int start_value){
    start = start_value;
};

void interval::set_end(int end_value){
    end = end_value;
};

void interval::set_interval(int start_value, int end_value){
    start = start_value;
    end = end_value;
};

int interval::get_start(void){
    return(start);
};

int interval::get_end(void){
    return(end);
};

bool interval::range_correct(void){
    if(interval::get_start() > interval::get_end()){
        return(true);
    };
    return(false);
};

int interval::check_overlap(interval comparator, bool only_overlap){ 
    /*
    1: overlaps,
    2: interval equal
    3:interval left_flank, 4: interval left_flank_overlap,    
    5: interval right_flank, 6: interval right_flank_overlap, 
    7: interval within left, 8: itnerval within right
    9: comparator within left, 10: comparator within right
    */
    const int comparator_start = comparator.get_start();
    const int comparator_end = comparator.get_end();
    const int interval_start = interval::get_start();
    const int interval_end = interval::get_end();

    bool start_leading = interval_start < comparator_start;
    bool start_lagging = interval_start > comparator_start;
    bool start_equal = interval_start == comparator_start;

    bool end_leading = interval_end < comparator_end;
    bool end_lagging = interval_end > comparator_end;
    bool end_equal = interval_end == comparator_end;

    bool start_end_leading = interval_start < comparator_end;
    bool start_end_lagging = interval_start > comparator_end;
    //bool start_end_equal = interval_start == comparator_end;

    bool end_start_leading = interval_end < comparator_start;
    bool end_start_lagging = interval_end > comparator_start;
   // bool end_start_equal = interval_end == comparator_start;

    bool interval_left_flank_overlap = start_leading && end_leading && end_start_lagging;
    bool interval_left_flank = start_leading && end_start_leading;
    bool interval_right_flank = start_end_lagging;
    bool interval_right_flank_overlap = start_lagging && end_lagging && start_end_leading;
    bool interval_equal = start_equal && end_equal;
    bool interval_within_left = start_equal && end_start_leading;
    bool interval_within_right = start_lagging && end_equal;
    bool comparator_within_left = start_equal && end_lagging;
    bool comparator_within_right = start_leading && end_equal;
    bool quit_after_overlap=only_overlap || interval_left_flank_overlap || interval_right_flank_overlap || interval_equal || interval_within_left || interval_within_right || comparator_within_left || comparator_within_right;

    if (quit_after_overlap){
        return(1);
    } else {
        if (interval_equal){
            return(2);
        } else if (interval_left_flank){
            return(3);
        } else if (interval_left_flank_overlap){
            return(4);
        } else if (interval_right_flank){
            return(5);
        } else if (interval_right_flank_overlap){
            return(6);
        } else if (interval_within_left){
            return(7);
        } else if (interval_within_right){
            return(8);
        } else if (comparator_within_left){
            return(9);
        } else if (comparator_within_right){
            return(10);
        } else {
            return(0);
        };
    };
};

bool interval::operator== (interval interval_object){
    bool start = interval::get_start() == interval_object.get_start();
    bool end = interval::get_end()==interval_object.get_end();
    if (start && end){
        return(true);
    } else {
        return(false);
    }
};

/********************************************************************************************************/
// The ginterval class function definitions

ginterval::ginterval(){};

ginterval::ginterval(int start_value, int end_value, std::string strand_value, std::string chromosome_value){
    start = start_value;
    end = end_value;
    strand = strand_value;
    chromosome = chromosome_value;
};

void ginterval::set_strand(std::string strand_value){
    strand = strand_value;
};

std::string ginterval::get_strand(void){
    return(strand);
};

void ginterval::set_chromosome(std::string chromosome_value){
    chromosome = chromosome_value;
};

std::string ginterval::get_chromosome(void){
    return(chromosome);
};

bool compare_gintervals(ginterval first, ginterval second){ // must be on the same chromosome for this to work
    // https://stackoverflow.com/q/1380463

    if (first.get_chromosome() == second.get_chromosome()){
        if (first.get_strand()=="-"){
            if (first.get_start() > second.get_start()){
                return(true);
            } else {
                return(false);
            }
        } else {
            if (first.get_start() < second.get_start()){
                return(true);
            } else {
                return(false);
            }
        }
    } else {
        return(false);
    }
}

bool ginterval::operator== (ginterval ginterval_object){
    bool start = ginterval::get_start() == ginterval_object.get_start();
    bool end = ginterval::get_end()==ginterval_object.get_end();
    bool chromosome = ginterval::get_chromosome()==ginterval_object.get_chromosome();
    bool strand = ginterval::get_strand()==ginterval_object.get_strand();
    if (start && end && chromosome && strand){
        return(true);
    } else {
        return(false);
    }
};


/********************************************************************************************************/
// The gset class function definitions

gset::gset(void){};

gset::gset(std::vector<ginterval> ginterval_vector,std::string parent_value){
    gintervals = ginterval_vector;
    parent = parent_value;
};

void gset::sort_gintervals(void){ // sorts gintervals by chromosome, start, end; must be on same strand for this to work
    // https://stackoverflow.com/q/1380463

    sort(gintervals.begin(),gintervals.end(),compare_gintervals);
};

std::vector<ginterval> gset::get_gintervals(void){
    return(gintervals);
};

void gset::set_parent(std::string parent_value){
    parent = parent_value;
};

void gset::add_ginterval(ginterval ginterval_to_add){
    gintervals.push_back(ginterval_to_add);
};

ginterval gset::get_range(void){ // gset must have non-overlapping gintervals, same chromosome, same strand
    std::string strand = gintervals[0].get_strand();
    std::string chromosome = gintervals[0].get_chromosome();
    int start = gset::gintervals[0].get_start();
    int end = gset::gintervals[gset::gintervals.size()-1].get_end();
    return(ginterval(start,end,strand,chromosome));
};

/********************************************************************************************************/
// The transcript class function definitions

transcript::transcript(void){};

transcript::transcript(gset exons_gset, gset cds_gset, gset utr5_gset, gset utr3_gset){
    exons = exons_gset;
    cds = cds_gset;
    utr5 = utr5_gset;
    utr3 = utr3_gset;
};

bool transcript::add_exon(ginterval exon_to_add){ // returns false if fail, true if succeed in adding
    bool equality_constraint=true;
    std::vector<ginterval> exon_vector = exons.get_gintervals();
    for (ginterval exon : exon_vector){
        if (exon == exon_to_add){
            equality_constraint=false;
            break;
        };
    };
    // if there are any exons equal to the one being added, don't add it.
        if (equality_constraint==true){
            exons.add_ginterval(exon_to_add);
        }
        return(equality_constraint);
};

bool transcript::add_cds(ginterval cds_to_add){ // returns false if fail, true if succeed in adding
    bool equality_constraint=true;
    std::vector<ginterval> cds_vector = cds.get_gintervals();
    for (ginterval single_cds : cds_vector){
        if (single_cds == cds_to_add){
            equality_constraint=false;
            break;
        };
    };
    if (equality_constraint==true){
        cds.add_ginterval(cds_to_add);
    }
    return(equality_constraint);
    // if there are any exons equal to the one being added, don't add it.
};

bool transcript::add_utr5(ginterval utr5_to_add){
    bool equality_constraint=true;
    std::vector<ginterval> utr5_vector = utr5.get_gintervals();
    for (ginterval single_utr5 : utr5_vector){
        if (single_utr5 == utr5_to_add){
            equality_constraint=false;
            break;
        };
    };
    if (equality_constraint==true){
        utr5.add_ginterval(utr5_to_add);
    }
    return(equality_constraint);
};

bool transcript::add_utr3(ginterval utr3_to_add){
    bool equality_constraint=true;
    std::vector<ginterval> utr3_vector = utr3.get_gintervals();
    for (ginterval single_utr3 : utr3_vector){
        if (single_utr3 == utr3_to_add){
            equality_constraint=false;
            break;
        };
    };
    if (equality_constraint==true){
        utr5.add_ginterval(utr3_to_add);
    }
    return(equality_constraint);
};

gset transcript::get_exons(void){
    return(exons);
};

gset transcript::get_cds(void){
    return(cds);
};

gset transcript::get_utr5(void){
    return(utr5);
};

gset transcript::get_utr3(void){
    return(utr3);
};

std::string transcript::get_identity(void){
    return(identity);
};

bool is_protein_coding(void){return(true);};

bool is_valid(void){return(true);};


/********************************************************************************************************/
// The gene class function definitions

gene::gene(void){};

gene::gene(int start_value,int end_value, std::string strand_value, std::string chromosome_value,
std::vector<transcript> transcripts_vector, std::string identity_value){
    start = start_value;
    end = end_value;
    strand = strand_value;
    chromosome = chromosome_value;
    transcripts = transcripts_vector;
    identity = identity_value;
};

int gene::get_start(void){
    return(start);
};

int gene::get_end(void){
    return(end);
};

std::string gene::get_strand(void){
    return(strand);
};

std::string gene::get_identity(void){
    return(identity);
}

std::string gene::get_chromomosome(void){
    return(chromosome);
};

void gene::add_transcript(transcript transcript_to_add){
    transcripts.push_back(transcript_to_add);
    transcript_identities.push_back(transcript_to_add.get_identity());
};

std::vector<std::string> gene::get_transcript_identities(void){
    return(transcript_identities);
};

int gene::find_transcript(std::string transcript_identity){
    for (size_t i=0; i<transcript_identities.size(); i++){
        if (transcript_identity == transcript_identities[i]){
            return(i);
        }
    }
    return(-1);
};

transcript gene::get_transcript(int tx_index){ // mean to be used after find_transcript with external error handling
    return(transcripts[tx_index]);
};

bool gene::is_emtpy(void){
    return(identity.empty());
};
