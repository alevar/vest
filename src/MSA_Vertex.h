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
    MSA_Vertex() = default;
    explicit MSA_Vertex(int num_ref, int pos){
        this->contents = std::vector<uint8_t>(num_ref);
        this->pos = pos;
    };
    ~MSA_Vertex() = default;

    void add_snp(uint8_t nt, uint16_t ref){
        this->contents[ref] = nt;
    }

    std::vector<uint8_t> get_contents() const {
        return this->contents;
    }

    int get_pos() const {
        return this->pos;
    }

    typedef std::vector<uint8_t>::iterator iterator;
    typedef std::vector<uint8_t>::const_iterator const_iterator;
    typedef std::vector<uint8_t>::reference reference;
    iterator begin() {return contents.begin();}
    const_iterator begin() const { return contents.begin();}
    iterator end() {return contents.end();}
    const_iterator end() const { return contents.end();}

    bool operator==(const MSA_Vertex& mv) const{
        return this->contents == mv.get_contents() && this->pos == mv.get_pos();
    }

private:
    std::vector<uint8_t> contents; // vector of nucleotides, where nucleotide is kept at the position of the reference id
    uint32_t pos = 0; // position of the vertex in the MSA

};

#endif //VEST_MSA_VERTEX_H
