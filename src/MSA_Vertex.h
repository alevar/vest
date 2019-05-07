//
// Created by sparrow on 5/7/19.
//

#ifndef VEST_MSA_VERTEX_H
#define VEST_MSA_VERTEX_H


#include <cstdint>
#include <vector>
#include <unordered_map>

class MSA_Vertex {
public:
    explicit MSA_Vertex(int num_ref){
        this->contents = std::vector<uint8_t>(num_ref);
    };
    ~MSA_Vertex() = default;

    void add_snp(uint8_t nt, uint16_t ref){
        this->contents[ref] = nt;
    }

    std::vector<uint8_t> get_contents() const {
        return this->contents;
    }

    typedef std::vector<uint8_t>::iterator iterator;
    typedef std::vector<uint8_t>::const_iterator const_iterator;
    typedef std::vector<uint8_t>::reference reference;
    iterator begin() {return contents.begin();}
    const_iterator begin() const { return contents.begin();}
    iterator end() {return contents.end();}
    const_iterator end() const { return contents.end();}

    bool operator==(const MSA_Vertex& mv) const{
        return this->contents == mv.get_contents();
    }

private:
    std::vector<uint8_t> contents; // vector of nucleotides, where nucleotide is kept at the position of the reference id

};

#endif //VEST_MSA_VERTEX_H
