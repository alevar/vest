//
// Created by sparrow on 5/7/19.
//

#ifndef VEST_MSA_VERTEX_H
#define VEST_MSA_VERTEX_H

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <set>
#include <map>

class MSA_Vertex {
public:
    MSA_Vertex(){
        this->empty = true;
    };
    explicit MSA_Vertex(int num_ref, int pos){
        this->empty = false;
        this->contents = std::vector<std::tuple<std::string,uint32_t,uint16_t>>(num_ref,std::make_tuple("",0,0));
        this->pos = pos;
    };
    ~MSA_Vertex() = default;

    void add_snp(std::string nt, uint16_t ref){
        std::get<0>(this->contents[ref]) += nt;
    }

    void add_edge(uint32_t next,uint16_t ref){
        std::get<1>(this->contents[ref]) = next;
    }

    std::vector<std::tuple<std::string,uint32_t,uint16_t>> get_contents() const {
        return this->contents;
    }

    int get_pos() const {
        return this->pos;
    }

    std::string get_nt(uint16_t ref_id) const {
        return std::get<0>(this->contents[ref_id]);
    }

    bool isEmpty(){
        return this->empty;
    }

    typedef std::vector<std::tuple<std::string,uint32_t,uint16_t> >::iterator iterator;
    typedef std::vector<std::tuple<std::string,uint32_t,uint16_t>>::const_iterator const_iterator;
    typedef std::vector<std::tuple<std::string,uint32_t,uint16_t>>::reference reference;
    iterator begin() {return contents.begin();}
    const_iterator begin() const { return contents.begin();}
    iterator end() {return contents.end();}
    const_iterator end() const { return contents.end();}

    bool operator==(const MSA_Vertex& mv) const{
        return this->contents == mv.get_contents() && this->pos == mv.get_pos();
    }

    void save(std::ofstream& out_fp){
        out_fp << this->pos << "\t";

        for(int i=0;i<this->contents.size();i++){
            if(!std::get<0>(this->contents[i]).empty()){ // only store those that are not empty
                out_fp << i << ":" << std::get<0>(this->contents[i])<<";";
                if(std::get<1>(this->contents[i]) == 0){
                    out_fp << "0";
                }
                else{
                    out_fp << std::get<1>(this->contents[i]);
                }
                if(i < this->contents.size()){
                    out_fp << "\t";
                }
            }
        }
        out_fp << std::endl;
    }

    void load(std::string& out_fp){

    }

    void get_nt_string(std::string& res){
        std::set<char> res_nts;
        for(auto& v : this->contents){
            if(!std::get<0>(v).empty()){
                for(auto &n : std::get<0>(v)) {
                    res_nts.insert(n);
                }
            }
        }
        for(auto& v : res_nts){
            res+=v;
        }
    }

    int get_most_abundant_refID(int& refID){
        int count=0;
        int refid=0;
        for(int i=0;i<this->contents.size();i++){
            if(count<std::get<2>(this->contents[i])){
                count=std::get<2>(this->contents[i]);
                refID=i;
            }
        }
        return count;
    }

    // returns true if the vertex describes as a nucleotide on the given reference
    bool has_ref(int refID){
        if(std::get<0>(this->contents[refID]).empty()){
            return false;
        }
        return true;
    }

    bool is_mapped(){return this->mapped;}

    bool set_mapped(){this->mapped=true;}

    void get_supported_nt_string(std::string& res){
        if(!is_mapped()){
            get_nt_string(res);
        }
        else {
            // get the most abundant references first
            int max_abund = 0;
            std::set<char> res_nts;
            for (auto &rc : this->contents) {
                if (std::get<2>(rc) > max_abund) {
                    max_abund = std::get<2>(rc);
                    res_nts.clear();
                    if (!std::get<0>(rc).empty()) {
                        for (auto &n : std::get<0>(rc)) {
                            res_nts.insert(n);
                        }
                    }
                    else{
                        std::cerr<<"reference sequence not found"<<std::endl;
                        exit(-1);
                    }
                }
            }
            for (auto &v : res_nts) {
                res += v;
            }
        }
    }

    void get_next_vts(std::vector<int>& next_vts){
        std::set<int> res_vts;
        for(auto& v : this->contents){
            if(!std::get<0>(v).empty()){
                res_vts.insert(std::get<1>(v));
            }
        }
        for(auto& v : res_vts){
            next_vts.push_back(v);
        }
    }

    int get_next_pos4ref(int refID){
        return std::get<1>(contents[refID]);
    }

    void inc_ref(int refID){
        std::get<2>(contents[refID])++;
    }

private:
    // lookup table where
    // A/a - 0
    // C/c - 1
    // G/g - 2
    // T/t - 3
    // other codes may be supported from here as well
    uint8_t nt_table[256]= {99, 99, 99, 99, 99, 99, 99, 99,//0
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

    std::vector<std::tuple<std::string,uint32_t,uint16_t> > contents = {}; // vector of nucleotides, where nucleotide is kept at the position of the reference id. each position also stores the index of the next vertex in case there exists an edge. The last element keeps count of the number of times the reference is confirmed
    uint32_t pos = 0; // position of the vertex in the MSA
    bool empty = true;
    bool mapped = false;
};

#endif //VEST_MSA_VERTEX_H
