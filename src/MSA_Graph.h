//
// Created by sparrow on 5/2/19.
//

#ifndef VEST_MSA_GRAPH_H
#define VEST_MSA_GRAPH_H

#include <string>

#include "MSA_Vertex.h"
#include "MSA_Edge.h"
#include "MSA_Index.h"
#include "MSA_List.h"


class MSA_Graph {
public:
    MSA_Graph() = default;
    MSA_Graph(int length, int num_refs);
    ~MSA_Graph() = default;

    uint16_t add_ref(std::string& ref_name);
    void add_pos(uint16_t id,uint32_t old_pos,uint32_t new_pos);
    void add_snp(std::string nt,uint32_t pos,uint16_t ref_id);
    void add_edge(uint32_t prev,uint32_t next, uint16_t ref_id);
    MSA_Vertex* get_vertex(uint32_t pos);

    std::string get_id(uint16_t id);
    int get_num_refs();
    int get_len();

    std::string get_nt(uint32_t vt_pos,uint16_t ref_id);
    void save_index(std::ofstream& out_fp);
    void save_graph(std::ofstream& out_fp); // saves vertices

private:
    MSA_Index index; // index which holds ref IDs
    int length = 0;
    int num_refs = 0;

    MSA_List<MSA_Vertex> vertices;
};


#endif //VEST_MSA_GRAPH_H
