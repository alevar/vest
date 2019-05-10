//
// Created by sparrow on 5/2/19.
//

#include "MSA.h"
#include "FastaTools.h"

MSA::MSA(std::string msa_fname) {
    this->msa_fname = msa_fname;
    this->parse_msa();
}

void MSA::parse_msa() {
    // now test that the file still exists
    std::ifstream msa_fp;
    msa_fp.open(this->msa_fname.c_str(),std::ios::in);
    if(!msa_fp.good()){
        std::cerr<<"Multiple Sequence Alignment fasta file not found"<<std::endl;
        exit(-1);
    }
    msa_fp.close();
    std::cerr << "Reading from the Multiple Sequence Alignment file: " << this->msa_fname << std::endl;

    FastaReader fastaReader(this->msa_fname);
    FastaRecord cur_genome;

    // first get basic info from the references
    bool first_genome = true;
    while (fastaReader.good()) {
        fastaReader.next(cur_genome);
        // first make sure that the contig info is added to the contig id map along with the length
        this->num_refs++;
//        std::cout<<cur_genome.id_<<std::endl;
        if (first_genome){
            first_genome = false;
            if (cur_genome.seq_.empty()){
                std::cerr<<"the input length is 0"<<std::endl;
                exit(1);
            }
            this->msa_len = cur_genome.seq_.size();
        }
        if (this->msa_len != cur_genome.seq_.size()){
            std::cerr<<"Genomes have different length"<<std::endl;
            exit(1);
        }
    }

    this->graph = MSA_Graph(this->msa_len,this->num_refs); // initiate the graph

//    std::cerr<<"==============================="<<std::endl;

    // now need to add sequence data to the vertices
    fastaReader.reset();
    cur_genome.clear();
    uint16_t cur_genome_id;
    std::string cur_nt,cur_iupac;
    while (fastaReader.good()) {
        fastaReader.next(cur_genome);

        cur_genome_id = this->graph.add_ref(cur_genome.id_); // add reference to the index

//        std::cout<<cur_genome.id_<<std::endl;
        MSA_Vertex* last_mv; // last vertex for connecting edges
        bool unassigned_last = true;
        uint32_t old_seq_idx = 0; // keeps track of the base position within the original reference sequence

        for (int i=0; i<cur_genome.seq_.size();i++){
            cur_nt = cur_genome.seq_[i];
            if(cur_nt != "-"){
                this->graph.add_pos(cur_genome_id,old_seq_idx,i); // add old and new positions to the index
                old_seq_idx++;
                for (auto cur_iupac_nt : this->IUPAC_REV[cur_nt]){
                    this->graph.add_snp(std::string(1,cur_iupac_nt),i,cur_genome_id);
//                    std::cerr<<cur_nt<<"\t"<<std::string(1,cur_iupac_nt)<<std::endl;
                }

                if (cur_genome.seq_[i-1] == '-' && !unassigned_last){ // add edge between last vertex and new one
                    last_mv->add_edge(i,cur_genome_id);
                }
                else{
                    if (cur_genome.seq_[i-1] != '-' && i != 0){
                        this->graph.add_edge(i-1,1,cur_genome_id);
                    }
                }
            }
            else{
                if (cur_genome.seq_[i-1] != '-' && i != 0){
                    last_mv = this->graph.get_vertex(i-1);
                    unassigned_last = false;
                }
            }
        }
    }
}

// this function outputs the MSA into a file
void MSA::to_msa(std::string out_msa_fname) {
    std::ofstream msa_fp(out_msa_fname.c_str());

    int glen = this->graph.get_len();
    int gnum_ref = this->graph.get_num_refs();

    MSA_Vertex cur_vt;

    for(int n=0;n<gnum_ref;n++){
        msa_fp << ">" <<this->graph.get_id(n) << std::endl; // save name

        for(int i=0;i<glen;i++){ // output the sequence of the current reference genome
            if (i%60 == 0 && i != 0){ // break sequence at certain number of characters
                msa_fp<<std::endl;
            }
            this->IUPAC_it = this->IUPAC.find(this->graph.get_nt(i,n));
            if(this->IUPAC_it != this->IUPAC.end()){
                msa_fp<<this->IUPAC_it->second;
            }
            else{
                msa_fp<<"-";
            }

        }
        msa_fp<<std::endl;
    }

    msa_fp.close();
}

// this function outputs the FASTA sequences into a file without MSA information
void MSA::to_fasta(std::string out_msa_fname) {
    std::ofstream msa_fp(out_msa_fname.c_str());

    int glen = this->graph.get_len();
    int gnum_ref = this->graph.get_num_refs();

    MSA_Vertex cur_vt;

    for(int n=0;n<gnum_ref;n++){
        msa_fp << ">" <<this->graph.get_id(n) << std::endl; // save name

        int line_count = 0;
        bool line_started = false;
        for(int i=0;i<glen;i++){ // output the sequence of the current reference genome
            if (line_count%60 == 0 && line_started){ // break sequence at certain number of characters
                msa_fp<<std::endl;
                line_started = false;
            }
            this->IUPAC_it = this->IUPAC.find(this->graph.get_nt(i,n));
            if(this->IUPAC_it != this->IUPAC.end()){
                msa_fp<<this->IUPAC_it->second;
                line_count++;
                line_started = true;
            }
        }
        msa_fp<<std::endl;
    }

    msa_fp.close();
}

void MSA::save_graph(std::string out_graph_base_name) {

    if(out_graph_base_name.rfind('/')==out_graph_base_name.length()-1){
        out_graph_base_name.pop_back();
    }

    MSA::_save_graph(out_graph_base_name);
    MSA::save_graph_info(out_graph_base_name);
    MSA::save_graph_contig_info(out_graph_base_name);

    std::string msa_fname(out_graph_base_name);
    msa_fname.append("/db.mus");
    MSA::to_msa(msa_fname);

    std::string fasta_fname(out_graph_base_name);
    fasta_fname.append("/db.fasta");
    MSA::to_fasta(fasta_fname);
}


void MSA::serialize() {

}

// save the actual graph datastructure
void MSA::_save_graph(std::string out_base){
    std::string graph_fname(out_base);
    graph_fname.append("/db.graph");
    std::ofstream graph_fp(graph_fname.c_str());

    this->graph.save_graph(graph_fp);

    graph_fp.close();
}

// save general info such as length, number of references, etc
void MSA::save_graph_info(std::string out_base){

    std::string graph_fname(out_base);
    graph_fname.append("/db.info");
    std::ofstream graph_fp(graph_fname.c_str());

    graph_fp << this->graph.get_len() << std::endl;
    graph_fp << this->graph.get_num_refs() << std::endl;

    graph_fp.close();
}

// save information with contig names and corresponding ids
void MSA::save_graph_contig_info(std::string out_base){
    std::string graph_fname(out_base);
    graph_fname.append("/db.contig");
    std::ofstream graph_fp(graph_fname.c_str());

    // saves the MSA_Index from MSA_Graph
    this->graph.save_index(graph_fp);

    graph_fp.close();
}