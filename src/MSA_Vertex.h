//
// Created by sparrow on 5/7/19.
//

#ifndef VEST_MSA_VERTEX_H
#define VEST_MSA_VERTEX_H

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_map>

class MSA_Vertex {
public:
    MSA_Vertex(){
        this->empty = true;
    };
    explicit MSA_Vertex(int num_ref, int pos){
        this->empty = false;
        this->contents = std::vector<std::pair<std::string,uint32_t>>(num_ref);
        this->pos = pos;
    };
    ~MSA_Vertex() = default;

    void add_snp(std::string nt, uint16_t ref){
        if (this->contents[ref].first.empty()){
            this->contents[ref].first = "";
        }
        this->contents[ref].first += nt;
//        std::cerr<<this->pos<<"\t"<<ref<<"\t"<<this->contents[ref].first<<"\t"<<nt<<std::endl;
    }

    void add_edge(uint32_t next,uint16_t ref){
        this->contents[ref].second = next;
    }

    std::vector<std::pair<std::string,uint32_t>> get_contents() const {
        return this->contents;
    }

    int get_pos() const {
        return this->pos;
    }

    std::string get_nt(uint16_t ref_id) const {
//        std::cerr<<ref_id<<"\t"<<this->contents[ref_id].first<<std::endl;
        return this->contents[ref_id].first;
    }

    bool isEmpty(){
        return this->empty;
    }

    typedef std::vector<std::pair<std::string,uint32_t> >::iterator iterator;
    typedef std::vector<std::pair<std::string,uint32_t>>::const_iterator const_iterator;
    typedef std::vector<std::pair<std::string,uint32_t>>::reference reference;
    iterator begin() {return contents.begin();}
    const_iterator begin() const { return contents.begin();}
    iterator end() {return contents.end();}
    const_iterator end() const { return contents.end();}

    bool operator==(const MSA_Vertex& mv) const{
        return this->contents == mv.get_contents() && this->pos == mv.get_pos();
    }

private:
    // lookup table where
    // A/a - 0
    // C/c - 1
    // G/g - 2
    // T/t - 3
    // other codes may be supported from here as well
    uint8_t nt_table[256]={99, 99, 99, 99, 99, 99, 99, 99,//0
                            99, 99, 99, 99, 99, 99, 99, 99,//8
                            99, 99, 99, 99, 99, 99, 99, 99,//16
                            99, 99, 99, 99, 99, 99, 99, 99,//24
                            99, 99, 99, 99, 99, 99, 99, 99,//32
                            99, 99, 99, 99, 99, 99, 99, 99,//40
                            99, 99, 99, 99, 99, 99, 99, 99,//48
                            99, 99, 99, 99, 99, 99, 99, 99,//56
                            99, 0 , 99, 1 , 99, 99, 99, 2 ,//64
                            99, 99, 99, 99, 99, 99, 99, 99,//72
                            99, 99, 99, 99, 3 , 99, 99, 99,//80
                            99, 99, 99, 99, 99, 99, 99, 99,//88
                            99, 0 , 99, 1 , 99, 99, 99, 2 ,//96
                            99, 99, 99, 99, 99, 99, 99, 99,//104
                            99, 99, 99, 99, 3 , 99, 99, 99,//112
                            99, 99, 99, 99, 99, 99, 99, 99,//120
                            99, 99, 99, 99, 99, 99, 99, 99,//128
                            99, 99, 99, 99, 99, 99, 99, 99,//136
                            99, 99, 99, 99, 99, 99, 99, 99,//144
                            99, 99, 99, 99, 99, 99, 99, 99,//152
                            99, 99, 99, 99, 99, 99, 99, 99,//160
                            99, 99, 99, 99, 99, 99, 99, 99,//168
                            99, 99, 99, 99, 99, 99, 99, 99,//176
                            99, 99, 99, 99, 99, 99, 99, 99,//184
                            99, 99, 99, 99, 99, 99, 99, 99,//192
                            99, 99, 99, 99, 99, 99, 99, 99,//200
                            99, 99, 99, 99, 99, 99, 99, 99,//208
                            99, 99, 99, 99, 99, 99, 99, 99,//216
                            99, 99, 99, 99, 99, 99, 99, 99,//224
                            99, 99, 99, 99, 99, 99, 99, 99,//232
                            99, 99, 99, 99, 99, 99, 99, 99,//240
                            99, 99, 99, 99, 99, 99, 99, 99};//248

    // keep track of the number of times a given nucleotide is used
    // counters for each nucleotide should be incremented for each read that contains a given vertex
    // to be used during the realignment
    uint8_t nts[4]={0,0,0,0};

    std::vector<std::pair<std::string,uint32_t> > contents = {}; // vector of nucleotides, where nucleotide is kept at the position of the reference id. each position also stores the index of the next vertex in case there exists an edge
    uint32_t pos = 0; // position of the vertex in the MSA
    bool empty = true;

};

#endif //VEST_MSA_VERTEX_H
