//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_GRAPH_H
#define VEST_MSA_GRAPH_H

#include <cstddef>
#include <string>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "MSA_Vertex.h"
#include "MSA_Edge.h"
#include "MSA_Index.h"
#include "MSA_List.h"

#include "gff.h"
#include "GFaSeqGet.h"

#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

class MSA_Graph {
public:
    MSA_Graph();
    MSA_Graph(int length, int num_refs);
    ~MSA_Graph() = default;

    uint16_t add_ref(std::string ref_name);
    void add_ref(std::string ref_name,int ref_id);
    void add_pos(uint16_t id,uint32_t new_pos);
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
    void save_base_refs(std::string&out_fp); // for each mapped base - save the reference supporting that base

    void fit_read(int refID,int ref_start,int end,int& newStart, int& s, std::vector<int>& not_removed, std::unordered_set<int>& added,uint8_t* seq,int32_t qlen);
    void find_location(int refID, int ref_start, int end, int& new_start, int& s);

    int pre_fit_annotation(std::string in_gff);
    int fit_annotation(std::string in_gff,std::string out_gff);

    int get_gff_pos(int refID,int pos);

    int get_first_pos(int refID); // retrieves the first position in the graph where the reference begins
    int get_last_pos(int refID);
    void set_removed_cleanup(int start, int end);
    void set_removed_cleanup(int pos);
    void get_first_mapped_pos(int& pos,int& refID);
    void get_last_mapped_pos(int& pos,int&refID);
    int get_most_abundant_refID(int pos,int&refID);
    void add2refcount(int pos,int refID);
    void init_refcouts();

    void set_used(int refid);
    void set_used(int start, int end, int refid);

    int get_num_clean_removed(int pos);
    void clean_gaps(int start,int end); // first mapped base and last mapped base

private:
    MSA_Index index; // index which holds ref IDs
    int length = 0;
    int num_refs = 0;

    MSA_List<MSA_Vertex> vertices;

    // related to updating the vector of removed vertices
    std::vector<int> removed; // TODO: simple vector for now - needs to be replaced by a more robust solution
    std::vector<int> removed_cleanup; // these are the vertices that are removed during the cleanup

    int farthestEnd=0; // the value is set to the last position that was processed so far.
                        // Since reads are sorted with respect to the MSA, any nades prior to this value
                        // are not to be removed. Instead if a read demands a removal - the cigar of the read
                        // is to be modified with a respective insertion
    int memo_refID,memo_end;

    // IUPAC definitions
    std::unordered_map<std::string,std::string> IUPAC;

    std::unordered_map<std::string,std::string>::iterator IUPAC_it;

    std::unordered_map<std::string,std::string> IUPAC_REV;
};


#endif //VEST_MSA_GRAPH_H
