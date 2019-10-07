//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_GRAPH_H
#define VEST_MSA_GRAPH_H

#include <string>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "MSA_Vertex.h"
#include "MSA_Edge.h"
#include "MSA_Index.h"
#include "MSA_List.h"

#include "gff.h"
#include "GFaSeqGet.h"

class MSA_Graph {
public:
    MSA_Graph() = default;
    MSA_Graph(int length, int num_refs);
    ~MSA_Graph() = default;

    uint16_t add_ref(std::string ref_name);
    void add_ref(std::string ref_name,int ref_id);
    void add_pos(uint16_t id,uint32_t old_pos,uint32_t new_pos);
    void add_snp(std::string nt,uint32_t pos,uint16_t ref_id);
    void add_edge(uint32_t prev,uint32_t next, uint16_t ref_id);
    void add_vertex(int pos, MSA_Vertex mv);

    MSA_Vertex* get_vertex(uint32_t pos);

    std::string get_id(uint16_t id);
    int get_id(std::string id);
    int get_num_refs();
    int get_len();
    int get_new_position(std::string& ref_name,int pos);
    int get_new_position(int refID, int pos);

    std::string get_nt(uint32_t vt_pos,uint16_t ref_id);
    void save_index(std::ofstream& out_fp);
    void save_graph(std::ofstream& out_fp); // saves vertices
    void save_graph2dot(std::ofstream &out_fp); // save vertices as a dot file
    void save_merged_fasta(std::string& out_fp);

    void fit_read(int refID,int ref_start,int end,int& newStart, int& s, std::vector<int>& not_removed, std::vector<int>& added);
    void find_location(int refID, int ref_start, int end, int& new_start, int& s);

    void fit_annotation(std::string in_gff,std::string out_gff);

    int get_gff_pos(int refID,int pos);

private:
    MSA_Index index; // index which holds ref IDs
    int length = 0;
    int num_refs = 0;

    MSA_List<MSA_Vertex> vertices;

    // related to updating the vector of removed vertices
    std::vector<int> removed; // TODO: simple vector for now - needs to be replaced by a more efficient solution
    int farthestEnd=0; // the value is set to the last position that was processed so far.
                        // Since reads are sorted with respect to the MSA, any nades prior to this value
                        // are not to be removed. Instead if a read demands a removal - the cigar of the read
                        // is to be modified with a respective insertion
    int memo_refID,memo_end;

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


#endif //VEST_MSA_GRAPH_H
