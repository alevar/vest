//
// Created by sparrow on 5/2/19.
//

#include "MSA_Graph.h"

MSA_Graph::MSA_Graph(int length,int num_refs) {
    this->length = length;
    this->num_refs = num_refs;

    // initialize vertices
    for(int i=0;i<length;i++){
        vertices.insert(MSA_Vertex(num_refs,i));
    }
}

// this function adds a reference name to the index and create a unique ID
uint16_t MSA_Graph::add_ref(std::string& ref_name) {
    return this->index.addRef(ref_name);
}

// this function sets a mapping between old and new positions in the refence/MSA
void MSA_Graph::add_pos(uint16_t id, uint32_t old_pos, uint32_t new_pos) {
    this->index.add(id,old_pos,new_pos);
}

// this function sets a snp for a given vertex
void MSA_Graph::add_snp(std::string nt, uint32_t pos, uint16_t ref_id) {
    MSA_Vertex* mv = this->vertices.get(pos);
    mv->add_snp(nt, ref_id);
    mv = this->vertices.get(pos);
//    std::cerr<<"g\t"<<mv->get_nt(ref_id)<<std::endl;
//    std::cerr<<nt<<std::endl;
}

void MSA_Graph::add_edge(uint32_t prev, uint32_t next, uint16_t ref_id) {
    MSA_Vertex* mv = this->vertices.get(prev);
    mv->add_edge(next, ref_id);
}

MSA_Vertex* MSA_Graph::get_vertex(uint32_t pos) {
    return this->vertices.get(pos);
}

std::string MSA_Graph::get_id(uint16_t id) {
    return this->index.getRef(id);
}

int MSA_Graph::get_num_refs(){
    return this->num_refs;
}

int MSA_Graph::get_len(){
    return this->length;
}

std::string MSA_Graph::get_nt(uint32_t vt_pos,uint16_t ref_id){
//    std::cerr<<vt_pos<<"\t";
    MSA_Vertex* mv = this->vertices.get(vt_pos);
    return mv->get_nt(ref_id);
}