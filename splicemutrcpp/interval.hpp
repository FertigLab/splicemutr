/* 
a library for handling genomic intervals and sets of genomic intervals (ie transcripts and genes)
Created: 20220824
Theron Palmer Jr.
*/

#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <string>
#include <vector>
#include <algorithm>

class interval{
    protected:
        int start;
        int end;
    public:
        interval(void);
        interval(int start_value,int end_value);
        void set_start(int start_value);
        void set_end(int end_value);
        void set_interval(int start_value, int end_value);
        int get_start(void);
        int get_end(void);
        bool range_correct(void);
        int check_overlap(interval comparator, bool only_overlap); // only_overlap==true for quick check of overlap
        bool operator== (interval interval_object);
};

class ginterval : public interval{
    protected:
        std::string chromosome;
        std::string strand; // 1:+; 2:-; 3:*
    public:
        ginterval(void);
        ginterval(int start_value,int end_value, std::string strand_value,std::string chromosome);
        void set_strand(std::string strand_value);
        std::string get_strand(void);
        void set_chromosome(std::string chromosome_value);
        std::string get_chromosome(void);
        bool operator== (ginterval ginterval_object);

};

bool compare_gintervals(ginterval first, ginterval second);

class gset{
    protected:
        std::string parent;
        std::vector<ginterval> gintervals;
    public:
        gset(void);
        gset(std::vector<ginterval> ginterval_vector,std::string parent_value);
        void sort_gintervals(void); // sorts gintervals by chromosome, start, end, strand
        std::vector<ginterval> get_gintervals(void);
        void set_parent(std::string parent_value);
        void add_ginterval(ginterval ginterval_to_add);
        ginterval get_range(void);
};

class transcript{
    protected:
        std::string identity;
        gset exons;
        gset cds;
        gset utr5;
        gset utr3;
    public:
        transcript(void);
        transcript(gset exons_gset, gset cds_gset, gset utr5_gset, gset utr3_gset);
        bool add_exon(ginterval exon_to_add);
        bool add_cds(ginterval cds_to_add);
        bool add_utr5(ginterval utr5_to_add);
        bool add_utr3(ginterval utr3_to_add);
        gset get_exons(void);
        gset get_cds(void);
        gset get_utr5(void);
        gset get_utr3(void);
        std::string get_identity(void);
        bool is_protein_coding(void); // this needs work, need lookup table for codon to amino-acid
        bool is_valid(void); // what makes a transcript valid?
};

class gene{
    protected:
        std::string identity;
        int start;
        int end;
        std::string strand;
        std::string chromosome;
        std::vector<transcript> transcripts;
        std::vector<std::string> transcript_identities;
    public:
        gene(void);
        gene(int start_value,int end_value, std::string strand_value, std::string chromosome_value,
        std::vector<transcript> transcripts_vector, std::string identity_value);
        int get_start(void);
        int get_end(void);
        std::string get_strand(void);
        std::string get_identity(void);
        std::string get_chromomosome(void);
        void add_transcript(transcript transcript_to_add);
        std::vector<std::string> get_transcript_identities(void);
        int find_transcript(std::string transcript_identity);
        transcript get_transcript(int tx_index);
        bool is_emtpy(void);
};

#endif