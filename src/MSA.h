//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_H
#define VEST_MSA_H

#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <map>

#include <htslib/sam.h>

#include "MSA_Graph.h"
#include "MSA_Vertex.h"

#define MAX_CIGARS 1024

class MSA {
public:
    MSA() = default;
    MSA(std::string msa_fname);
    ~MSA() = default;

    void to_msa(std::string out_msa_fname);
    void to_fasta(std::string out_fasta_fname);

    void save_graph(std::string out_graph_fname,std::string cmd);

    void load_graph(std::string in_graph_fname, std::string cmd);

    void realign(std::string in_sam,std::string out_sam);

private:
    std::string msa_fname;
    std::string msa_header_fname;
    FILE* msa_fhandle;
    int msa_len = 0; // number of nucleotides in the msa
    int num_refs = 0;

    std::string build_cmd = ""; // command used for building the graph
    std::string realign_cmd = ""; // command used for realigning the graph

    MSA_Graph graph;

    bool isMod(bam1_t* in_rec);
    void split_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header);
    void write_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header);
    void add_orig_ref_tags(bam1_t* in_rec,int ref);
    void add_split_tags(bam1_t* in_rec,int cur_slice,int opcode);

    void clean(int first_pos);
    void parse_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header);
    void change_cigar(bam1_t* in_rec,int s);
    void create_del(bam1_t* in_rec,std::vector<int>& not_removed);
    void create_ins(bam1_t* in_rec,std::vector<int>& added);

    void parse_msa();
    void save_graph_info(std::string out_base);
    void save_graph_contig_info(std::string out_base);
    void _save_graph(std::string out_base);
    void generate_bam_header(std::string out_base);

    void load_graph_info(std::ifstream& stream);
    void load_graph_contig_info(std::ifstream& stream);
    void _load_graph(std::ifstream& stream);

    bool change_data(bam1_t *in_rec,int num_cigars,int* cigars,int cur_start,int cur_len,bool shift);


    // IUPAC definitions
    std::unordered_map<std::string,std::string> IUPAC = std::unordered_map<std::string,std::string>({{"A","A"},
                                                                                                     {"C","C"},
                                                                                                     {"G","G"},
                                                                                                     {"T","T"},
                                                                                                     {"N","N"},
                                                                                                     {"AG","R"},{"GA","R"},
                                                                                                     {"CT","Y"},{"TC","Y"},
                                                                                                     {"CG","S"},{"GC","S"},
                                                                                                     {"AT","W"},{"TA","W"},
                                                                                                     {"GT","K"},{"TG","K"},
                                                                                                     {"AC","M"},{"CA","M"},
                                                                                                     {"CGT","B"},{"CTG","B"},{"GTC","B"},{"GCT","B"},{"TCG","B"},{"TGC","B"},
                                                                                                     {"AGT","D"},{"ATG","D"},{"GAT","D"},{"GTA","D"},{"TAG","D"},{"TGA","D"},
                                                                                                     {"ACT","H"},{"ATC","H"},{"CAT","H"},{"CTA","H"},{"TCA","H"},{"TAC","H"},
                                                                                                     {"ACG","V"},{"AGC","V"},{"CAG","V"},{"CGA","V"},{"GCA","V"},{"GAC","V"},
                                                                                                     {"ACGT","N"},{"ACTG","N"},{"AGCT","N"},{"AGTC","N"},{"ATCG","N"},{"ATGC","N"},
                                                                                                     {"CAGT","N"},{"CATG","N"},{"CGAT","N"},{"CGTA","N"},{"CTAG","N"},{"CTGA","N"},
                                                                                                     {"GACT","N"},{"GATC","N"},{"GCAT","N"},{"GCTA","N"},{"GTAC","N"},{"GTCA","N"},
                                                                                                     {"TACG","N"},{"TAGC","N"},{"TCAG","N"},{"TCGA","N"},{"TGAC","N"},{"TGCA","N"}});

    std::unordered_map<std::string,std::string>::iterator IUPAC_it;

    std::unordered_map<std::string,std::string> IUPAC_REV = std::unordered_map<std::string,std::string>({{"A","A"},{"a","A"},
                                                                                                         {"C","C"},{"c","C"},
                                                                                                         {"G","G"},{"g","G"},
                                                                                                         {"T","T"},{"t","T"},
                                                                                                         {"R","AG"},{"r","AG"},
                                                                                                         {"Y","CT"},{"y","CT"},
                                                                                                         {"S","CG"},{"s","CG"},
                                                                                                         {"W","AT"},{"w","AT"},
                                                                                                         {"K","GT"},{"k","GT"},
                                                                                                         {"M","AC"},{"m","AC"},
                                                                                                         {"B","CGT"},{"b","CGT"},
                                                                                                         {"D","AGT"},{"d","AGT"},
                                                                                                         {"H","ACT"},{"h","ACT"},
                                                                                                         {"V","ACG"},{"v","ACG"},
                                                                                                         {"N","ACGT"},{"n","ACTG"}});

};


#endif //VEST_MSA_H
