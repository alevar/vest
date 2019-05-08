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

class MSA {
public:
    MSA() = default;
    MSA(std::string msa_fname);
    ~MSA() = default;

private:
    std::string msa_fname;
    FILE* msa_fhandle;
    int msa_len = 0; // number of nucleotides in the msa
    int num_refs = 0;

    MSA_Graph graph;

    void parse_msa();

};


#endif //VEST_MSA_H
