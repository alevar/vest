//
// Created by sparrow on 5/2/19.
//

#include "MSA.h"

MSA::MSA(std::string msa_fname) {
    this->msa_fname = msa_fname;
}

void MSA::parse_msa() {
    this->msa_fhandle = fopen(this->msa_fname.c_str(), "r");
    if (msa_fhandle == nullptr){
        std::cerr << "FATAL: Couldn't open Multiple Sequence Alignment file: " << this->msa_fhandle<< std::endl;
        exit(1);
    }
    std::cout << "Reading from the Multiple Sequence Alignment file: " << this->msa_fhandle << std::endl;


}
