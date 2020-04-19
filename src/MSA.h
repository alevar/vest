//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_H
#define VEST_MSA_H

#include <cstddef>
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
#include <unordered_set>
#include <unordered_map>

#include <htslib/sam.h>

#include "MSA_Graph.h"
#include "MSA_Vertex.h"

#define MAX_CIGARS 1024

class MSA {
public:
    MSA();
    MSA(std::string msa_fname);
    ~MSA() = default;

    void to_msa(std::string out_msa_fname);
    void to_fasta(std::string out_fasta_fname);

    void save_graph(std::string out_graph_fname,std::string cmd);

    void load_graph(std::string in_graph_fname, std::string cmd);

    void realign(std::string in_sam);

    void pre_fit_annotation(std::string in_gff);
    void fit_annotation(std::string in_gff, std::string out_gff);
    void fit_bed(std::string in_bed,std::string out_bed);

    void set_gapfillname(std::string ref_name);
    void set_projection_name(std::string ref_name);

    void set_tmp_dir(std::string tmp_dir);
    void remove_tmp();
    void init_tmp();
    void set_out_fname(std::string out_fname);

private:
    std::string out_fname;
    std::string base_name;
    std::string tmp_dir = "";
    bool keep_tmp_data = false;

    int gap_fillID = -1;

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
    void add_orig_ref_tags(bam1_t* in_rec,int ref,int new_end);
    void add_split_tags(bam1_t* in_rec,int cur_slice,int opcode,int length);
    void remove_aux_tags(bam1_t* rec);
    int get_ins_len(bam1_t* rec);
    void get_ins_seq(bam1_t* rec,uint8_t* seq,int& seq_len);
    void get_ins_qual(bam1_t* rec,uint8_t* qual,int& qual_len);
    void get_split_op_tags(bam1_t* rec,int& opcode,int&oplen);
    void add_seq_slice_tag(bam1_t* in_rec,uint8_t* seq_slice,int seq_len);
    bool get_seq_slice(bam1_t* in_rec,bam1_t* out_rec,int cur_start,int cur_len,bool shift);

    void clean();
    void parse_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header);
    void change_cigar(bam1_t* in_rec,int s);
    void l2range(std::vector<int>& l,std::vector<std::pair<int,int>>& r);
    void add_cigar(bam1_t *curAl,int num_cigars,int* cigars);
    void create_del(bam1_t* in_rec,std::vector<std::pair<int,int>>& not_removed);
    void create_ins_old(bam1_t* in_rec,std::vector<std::pair<int,int>>& added);
    void create_ins(bam1_t* in_rec,std::unordered_set<int>& added);

    void parse_msa();
    void save_graph_info(std::string out_base);
    void save_graph_contig_info(std::string out_base);
    void _save_graph(std::string out_base);
    void generate_bam_header(std::string out_base);

    void load_graph_info(std::ifstream& stream);
    void load_graph_contig_info(std::ifstream& stream);
    void _load_graph(std::ifstream& stream);

    bool change_data(bam1_t *in_rec,int num_cigars,int* cigars,int cur_start,int cur_len,bool shift);

    void join_cigars_old(std::vector<bam1_t*>& reads,uint8_t *cigars,int& new_n_cigar_bytes);
    void get_dist2next(bam1_t* rec,bam1_t* next,int opcode,int oplen);
    void join_cigars(std::vector<bam1_t*>& reads,uint8_t *cigars,int& new_n_cigar_bytes,std::unordered_set<int>& added);
    void merge_seqs(uint8_t* data,uint8_t* seq,int seq_len,int& cur_mem_pos,bool orphan);
    void join_seqs(std::vector<bam1_t*>& reads,uint8_t* data,int& cur_mem_pos,int max_seq_len);
    void join_quals(std::vector<bam1_t*>& reads,uint8_t* data,int& cur_mem_pos,int max_seq_len);
    void join_reads(std::vector<bam1_t*>& reads,samFile *outSAM_joined,bam_hdr_t *outSAM_joined_header,int mate_start);


    // IUPAC definitions
    std::unordered_map<std::string,std::string> IUPAC;

    std::unordered_map<std::string,std::string>::iterator IUPAC_it;

    std::unordered_map<std::string,std::string> IUPAC_REV;
};

#endif //VEST_MSA_H
