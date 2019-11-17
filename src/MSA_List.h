//
// Created by Ales Varabyou on 5/6/19.
// Container for a graph (vertices and edges)
//

#ifndef VEST_MSA_LIST_H
#define VEST_MSA_LIST_H

#include <unordered_map>
#include <vector>
#include <iostream>

#include "MSA_Vertex.h"

inline void hash_combine(std::size_t& seed, std::size_t v){
    seed ^= v + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

struct MSA_Vertex_Contents_Hasher{ // based on the contents of the vertex
    std::size_t operator()(const MSA_Vertex& mv) const{
        std::size_t hash = 0;
        for (std::tuple<std::string,uint32_t,uint16_t> i : mv){
            for (int i2=0;i2<std::get<0>(i).size();i2++) {
                hash_combine(hash, (uint8_t)std::get<0>(i)[i2]);
            }
        }
        return hash;
    }
};

struct MSA_Vertex_Pos_Hasher{ // based on the position of the vertex
    std::size_t operator()(const MSA_Vertex& mv) const{
        return mv.get_pos();
    }
};

template <class OBJ>
class MSA_List {
public:
    MSA_List() = default;

    ~MSA_List() = default;

    void insert(OBJ ob) {
        vkeys.push_back(ob);
        hk_it = hkeys.insert(std::make_pair(ob,vkeys.size()-1));
        if (!hk_it.second){
            std::cerr<<"key already exists!"<<std::endl;
        }
    }

    void change(int pos,OBJ ob) {
        vkeys[pos] = ob;
        hkeys[ob] = pos;
    }

    void remove(OBJ ob) {
        uint32_t idx = hkeys[ob];
        vkeys[idx] = vkeys[vkeys.size()-1];
        hkeys[vkeys.back()] = idx;
        hkeys.erase(ob);
        vkeys.pop_back();
    }

    OBJ* get(int idx) {
        return &(*(vkeys.begin()+idx));
    };

    int size(){
        return vkeys.size();
    }

private:
    std::unordered_map<OBJ,uint32_t,MSA_Vertex_Pos_Hasher> hkeys;
    std::vector<OBJ> vkeys;

    std::pair<typename std::unordered_map<OBJ,uint32_t,MSA_Vertex_Pos_Hasher>::iterator,bool> hk_it;

};

#endif //VEST_MSA_LIST_H
