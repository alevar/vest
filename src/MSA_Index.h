//
// Created by sparrow on 5/6/19.
//

#ifndef VEST_MSA_INDEX_H
#define VEST_MSA_INDEX_H

#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

class MSA_Index {
public:
    MSA_Index() = default;
    ~MSA_Index() = default;

    uint16_t add_ref(std::string ref){
        ri_it = ref_to_id.insert(std::make_pair(ref,maxID));
        if (ri_it.second){ // successfully inserted new reference sequence
            ir_it = id_to_ref.insert(std::make_pair(maxID,ref));
            if(!ir_it.second){
                std::cerr<<"something is wrong with the index. reference did not exist, but the id exists"<<std::endl;
            }
            else{
                pos_idx.insert(std::make_pair(maxID,std::vector<uint32_t>()));
            }
            maxID++;
        }
        else{
            std::cerr<<"detected duplicate sequences"<<std::endl;
        }
        return ri_it.first->second;
    }

    void add_ref(std::string ref_name, int ref_id){
        ri_it = ref_to_id.insert(std::make_pair(ref_name,ref_id));
        if (!ri_it.second){
            std::cerr<<"duplicate reference IDs detected"<<std::endl;
            exit(1);
        }
        id_to_ref.insert(std::make_pair(ref_id,ref_name));
    }

    uint16_t getID(std::string& ref){
        return ref_to_id[ref];
    }
    std::string getRef(uint16_t id){
        return id_to_ref[id];
    }

    void add(std::string& ref,uint32_t old_pos,uint32_t new_pos){
        __add(getID(ref),old_pos,new_pos);
    }

    void add(uint16_t id,uint32_t old_pos,uint32_t new_pos){
        __add(id,old_pos,new_pos);
    }

    uint16_t getNewPos(uint16_t id, uint16_t old_pos){
        return pos_idx[id][old_pos];
    }

    void save(std::ofstream& out_fp){
        for (auto& ri : this->ref_to_id){
            out_fp << ri.first << "\t" << ri.second << std::endl;
        }
    }

private:
    uint16_t maxID = 0;

    std::unordered_map<uint16_t,std::vector<uint32_t> > pos_idx; // for each position in the reference lists position within MSA
    std::pair<std::unordered_map<uint16_t,std::vector<uint32_t> >::iterator,bool> pi_it;

    std::unordered_map<std::string,uint16_t> ref_to_id;
    std::unordered_map<uint16_t,std::string> id_to_ref;

    std::pair<std::unordered_map<std::string,uint16_t>::iterator,bool> ri_it;
    std::pair<std::unordered_map<uint16_t,std::string>::iterator,bool> ir_it;

    void __add(uint16_t id,uint32_t old_pos,uint32_t new_pos){
        pos_idx[id].push_back(new_pos);
    }

};


#endif //VEST_MSA_INDEX_H
