//
// Created by sparrow on 5/2/19.
//

#include "MSA_Graph.h"

MSA_Graph::MSA_Graph(int length,int num_refs) {
    this->length = length;

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
void MSA_Graph::add_snp(const std::string nt, uint32_t pos, uint16_t ref_id) {
    MSA_Vertex mv = this->vertices.get(pos);
    mv.add_snp((uint8_t) nt[0], ref_id);
}
