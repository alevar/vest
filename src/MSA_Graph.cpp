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

uint16_t MSA_Graph::add_ref(std::string& ref_name) {
    return this->index.addRef(ref_name);
}
