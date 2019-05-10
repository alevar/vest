//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_H
#define VEST_MSA_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <map>

#include "MSA_Graph.h"
#include "MSA_Vertex.h"

class MSA {
public:
    MSA() = default;
    MSA(std::string msa_fname);
    ~MSA() = default;

    void to_msa(std::string out_msa_fname);
    void to_fasta(std::string out_fasta_fname);

    void save_graph(std::string out_graph_fname);

private:
    std::string msa_fname;
    FILE* msa_fhandle;
    int msa_len = 0; // number of nucleotides in the msa
    int num_refs = 0;

    MSA_Graph graph;

    void parse_msa();
    void serialize();
    void save_graph_info(std::string out_graph_info_fname);
    void save_graph_contig_info(std::string out_graph_contig_info_fname);
    void _save_graph(std::string out_graph_fname);


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
