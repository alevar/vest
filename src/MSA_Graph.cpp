//
// Created by sparrow on 5/2/19.
//

#include "MSA_Graph.h"

MSA_Graph::MSA_Graph(int length,int num_refs) {
    this->length = length;
    this->num_refs = num_refs;

    // initialize vertices
    for(int i=0;i<length;i++){
        vertices.insert(MSA_Vertex(num_refs,i));
    }
}

// this function adds a reference name to the index and create a unique ID
uint16_t MSA_Graph::add_ref(std::string ref_name) {
    return this->index.add_ref(ref_name);
}

void MSA_Graph::add_ref(std::string ref_name, int ref_id){
    this->index.add_ref(ref_name,ref_id);
}

// this function sets a mapping between old and new positions in the refence/MSA
void MSA_Graph::add_pos(uint16_t id, uint32_t old_pos, uint32_t new_pos) {
    this->index.add(id,old_pos,new_pos);
}

// this function sets a snp for a given vertex
void MSA_Graph::add_snp(std::string nt, uint32_t pos, uint16_t ref_id) {
    MSA_Vertex* mv = this->vertices.get(pos);
    mv->add_snp(nt, ref_id);
}

void MSA_Graph::add_edge(uint32_t prev, uint32_t next, uint16_t ref_id) {
    MSA_Vertex* mv = this->vertices.get(prev);
    mv->add_edge(next, ref_id);
}

MSA_Vertex* MSA_Graph::get_vertex(uint32_t pos) {
    return this->vertices.get(pos);
}

std::string MSA_Graph::get_id(uint16_t id) {
    return this->index.getRef(id);
}

int MSA_Graph::get_id(std::string id){
    return this->index.getID(id);
}

int MSA_Graph::get_num_refs() {
    return this->num_refs;
}

int MSA_Graph::get_len() {
    return this->length;
}

std::string MSA_Graph::get_nt(uint32_t vt_pos,uint16_t ref_id) {
    MSA_Vertex* mv = this->vertices.get(vt_pos);
    return mv->get_nt(ref_id);
}

int MSA_Graph::get_new_position(std::string &ref_name, int pos) {
    return this->index.getNewPos(ref_name,pos);
}

int MSA_Graph::get_new_position(int refID, int pos) {
    return this->index.getNewPos(refID,pos);
}

void MSA_Graph::save_index(std::ofstream& out_fp) {
    this->index.save(out_fp);
}

void MSA_Graph::save_graph(std::ofstream &out_fp) {
    MSA_Vertex* mv;
    for(int i=0;i<this->vertices.size();i++){
        mv = this->vertices.get(i);
        mv->save(out_fp);
    }
}

void MSA_Graph::save_graph2dot(std::ofstream &out_fp){
    MSA_Vertex* mv;
    MSA_Vertex* next_vt;
    std::string nt_str = "";
    std::string next_nt_str = "";
    std::vector<int> next_vts;
    for(int i=0;i<this->vertices.size();i++){
        mv = this->vertices.get(i);
        mv->get_nt_string(nt_str);
        mv->get_next_vts(next_vts);
        for(auto v : next_vts){
            next_vt = this->vertices.get(v);
            next_vt->get_nt_string(next_nt_str);
            out_fp<<i<<":"<<nt_str<<" -> "<<v<<":"<<next_nt_str<<";"<<std::endl;
            next_nt_str.clear();
        }
        nt_str.clear();
        next_vts.clear();
    }
}

void MSA_Graph::save_merged_fasta(std::string& out_fp){
    std::ofstream merged_fp(out_fp.c_str());
    merged_fp<<">MSA"<<std::endl;

    MSA_Vertex* mv;
    std::string iupac_nt;
    for(int i=0;i<this->vertices.size();i++){
        std::string nt_str = "";
        mv = this->vertices.get(i);
        mv->get_nt_string(nt_str);
        iupac_nt = this->IUPAC[nt_str];
        merged_fp<<iupac_nt;
        nt_str.clear();
        if((i+1) % 60 == 0){
            merged_fp<<std::endl;
        }
    }
    merged_fp<<std::endl;
    merged_fp.close();

    std::string merged_fai_fname(out_fp);
    merged_fai_fname.append(".fai");
    std::ofstream merged_fai_fp(merged_fai_fname);

    merged_fai_fp<<"MSA"<<"\t"<<this->vertices.size()<<"\t"<<5<<"\t"<<60<<"\t"<<61<<std::endl;

    merged_fai_fp.close();
}

void MSA_Graph::add_vertex(int pos,MSA_Vertex mv) {
    this->vertices.change(pos,mv);
}
