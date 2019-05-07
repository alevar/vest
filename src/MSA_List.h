//
// Created by Ales Varabyou on 5/6/19.
// Container for a graph (vertices and edges)
//

#ifndef VEST_MSA_LIST_H
#define VEST_MSA_LIST_H

#include <unordered_map>
#include <vector>
#include <iostream>

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

    void remove(OBJ ob) {
        uint32_t idx = hkeys[ob];
        vkeys[idx] = vkeys[vkeys.size()-1];
        hkeys[vkeys.back()] = idx;
        hkeys.erase(ob);
        vkeys.pop_back();
    }

    OBJ get(int idx) {
        return vkeys[idx];
    };

    int size(){
        return vkeys.size();
    }

private:
    std::unordered_map<OBJ,uint32_t> hkeys;
    std::vector<OBJ> vkeys;

    std::pair<typename std::unordered_map<OBJ,uint32_t>::iterator,bool> hk_it;

};

#endif //VEST_MSA_LIST_H
