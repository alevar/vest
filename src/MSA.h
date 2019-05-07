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

class MSA {
public:
    MSA() = default;
    MSA(std::string msa_fname);
    ~MSA() = default;

    void parse_msa();

private:
    std::string msa_fname;
    FILE* msa_fhandle;
    int msa_len = 0; // number of nucleotides in the msa
    int num_refs = 0;

};


#endif //VEST_MSA_H
