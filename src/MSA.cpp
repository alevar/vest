//
// Created by sparrow on 5/2/19.
//

#include "MSA.h"
#include "FastaTools.h"

void print_cigar(bam1_t *al){
    for (uint8_t c=0;c<al->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(al);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
        std::cout<<length<<bam_cigar_opchr(opcode);
    }
    std::cout<<std::endl;
}

void print_raw_cigar(int cigars[MAX_CIGARS],int n_cigar){
    for (uint8_t c=0;c<n_cigar;++c){
        int opcode=bam_cigar_op(cigars[c]);
        int length=bam_cigar_oplen(cigars[c]);
        std::cout<<length<<bam_cigar_opchr(opcode);
    }
    std::cout<<std::endl;
}

void print_raw_cigar_uint(uint32_t cigars[MAX_CIGARS],int n_cigar){
    for (uint8_t c=0;c<n_cigar;++c){
        int opcode=bam_cigar_op(cigars[c]);
        int length=bam_cigar_oplen(cigars[c]);
        std::cout<<length<<bam_cigar_opchr(opcode);
    }
    std::cout<<std::endl;
}

void print_seq(bam1_t *new_rec){
    int32_t qlen = new_rec->core.l_qseq;
    int8_t *buf = NULL;
    buf = static_cast<int8_t *>(realloc(buf, qlen+1));
    buf[qlen] = '\0';
    uint8_t* seq = bam_get_seq(new_rec);
    for (int i = 0; i < qlen; ++i)
        buf[i] = bam_seqi(seq, i);
    for (int i = 0; i < qlen; ++i) {
        buf[i] = seq_nt16_str[buf[i]];
    }
    std::string str_seq((char*)(char*)buf);
    std::cout<<str_seq<<std::endl;
}

void print_qual(bam1_t *new_rec){
    std::cout<<bam_get_qual(new_rec)<<std::endl;
}

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
    // now need to add sequence data to the vertices
    fastaReader.reset();
    cur_genome.clear();
    uint16_t cur_genome_id;
    std::string cur_nt,cur_iupac;
    while (fastaReader.good()) {
        fastaReader.next(cur_genome);

        cur_genome_id = this->graph.add_ref(cur_genome.id_); // add reference to the index

        MSA_Vertex* last_mv; // last vertex for connecting edges
        bool unassigned_last = true;

        for (int i=0; i<cur_genome.seq_.size();i++){
            cur_nt = cur_genome.seq_[i];
            if(cur_nt != "-"){
                this->graph.add_pos(cur_genome_id,i); // add old and new positions to the index
                for (auto cur_iupac_nt : this->IUPAC_REV[cur_nt]){
                    this->graph.add_snp(std::string(1,cur_iupac_nt),i,cur_genome_id);
//                    std::cerr<<cur_nt<<"\t"<<std::string(1,cur_iupac_nt)<<std::endl;
                }

                if (cur_genome.seq_[i-1] == '-' && !unassigned_last){ // add edge between last vertex and new one
                    last_mv->add_edge(i,cur_genome_id);
                }
                else{
                    if (cur_genome.seq_[i-1] != '-' && i != 0){
                        this->graph.add_edge(i-1,i,cur_genome_id);
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

void MSA::save_graph(std::string out_graph_base_name, std::string cmd) {

    this->build_cmd = cmd;

    if(out_graph_base_name.rfind('/')==out_graph_base_name.length()-1){
        out_graph_base_name.pop_back();
    }

    MSA::_save_graph(out_graph_base_name);
    MSA::save_graph_info(out_graph_base_name);
    MSA::save_graph_contig_info(out_graph_base_name);
    MSA::generate_bam_header(out_graph_base_name);

    std::string msa_fname(out_graph_base_name);
    msa_fname.append("/db.mus");
    MSA::to_msa(msa_fname);

    std::string fasta_fname(out_graph_base_name);
    fasta_fname.append("/db.fasta");
    MSA::to_fasta(fasta_fname);

    // save merged reference
    std::string merged_fasta_fname(out_graph_base_name);
    merged_fasta_fname.append("/db.merged.fasta");
    this->graph.save_merged_fasta(merged_fasta_fname);
}

// save the actual graph datastructure
void MSA::_save_graph(std::string out_base){
    std::string graph_fname(out_base);
    graph_fname.append("/db.graph");
    std::ofstream graph_fp(graph_fname.c_str());

    std::string graph_dot_fname(out_base);
    graph_dot_fname.append("/db.graph.dot");
    std::ofstream graph_dot_fp(graph_dot_fname.c_str());

    this->graph.save_graph(graph_fp);
    this->graph.save_graph2dot(graph_dot_fp);

    graph_fp.close();
    graph_dot_fp.close();
}

// save general info such as length, number of references, etc
void MSA::save_graph_info(std::string out_base){

    std::string graph_fname(out_base);
    graph_fname.append("/db.info");
    std::ofstream graph_fp(graph_fname.c_str());

    graph_fp << this->graph.get_len() << std::endl;
    graph_fp << this->graph.get_num_refs() << std::endl;
    graph_fp << this->build_cmd << std::endl;

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

// generates a sample bam header for the realignment to the given MSA
void MSA::generate_bam_header(std::string out_base) {
    std::string bam_header_fname(out_base);
    bam_header_fname.append("/db.header.bam");
    std::ofstream bam_header_fp(bam_header_fname.c_str());

    bam_header_fp << "@HD\tVN:1.0\tSO:unsorted" << std::endl;
    bam_header_fp << "@SQ\tSN:MSA\tLN:" << this->msa_len << std::endl;
    bam_header_fp << "@PG\tID:Vest\tPN:Vest\tVN:1.0\tCL:" << this->build_cmd << std::endl;

    bam_header_fp.close();
}

void MSA::load_graph(std::string in_graph_fname, std::string cmd) {
    this->realign_cmd = cmd;

    // first check that everything is there
    // then load the components
    if(in_graph_fname.rfind('/')==in_graph_fname.length()-1){
        in_graph_fname.pop_back();
    }
    std::string graph_fname(in_graph_fname);
    graph_fname.append("/db.graph");
    std::string info_fname(in_graph_fname);
    info_fname.append("/db.info");
    std::string contig_info_fname(in_graph_fname);
    contig_info_fname.append("/db.contig");
    std::string mus_fname(in_graph_fname);
    mus_fname.append("/db.mus");
    std::string fasta_fname(in_graph_fname);
    fasta_fname.append("/db.fasta");
    this->msa_header_fname = std::string(in_graph_fname);
    this->msa_header_fname.append("/db.header.bam");

    std::ifstream info_fp;
    info_fp.open(info_fname.c_str(),std::ios::in);
    if(!info_fp.good()){
        std::cerr<<"FATAL: Couldn't open graph information data: "<<info_fname<<std::endl;
        exit(1);
    }

    std::cerr<<"@LOG:::: Begin loading graph"<<std::endl;
    MSA::load_graph_info(info_fp);
    info_fp.close();

    std::ifstream contig_info_fp;
    contig_info_fp.open(contig_info_fname.c_str(),std::ios::in);
    if(!contig_info_fp.good()){
        std::cerr<<"FATAL: Couldn't open contig information file: "<<contig_info_fname<<std::endl;
        exit(1);
    }
    MSA::load_graph_contig_info(contig_info_fp);
    contig_info_fp.close();

    std::ifstream graph_fp;
    graph_fp.open(graph_fname.c_str(),std::ios::in);
    if(!graph_fp.good()){
        std::cerr<<"FATAL: Couldn't open graph contents file: "<<graph_fname<<std::endl;
        exit(1);
    }
    MSA::_load_graph(graph_fp);
    graph_fp.close();

    // check the header file
    std::ifstream msa_header_fp;
    msa_header_fp.open(msa_header_fname.c_str(),std::ios::in);
    if(!msa_header_fp.good()){
        std::cerr<<"FATAL: Couldn't open msa header contents file: "<<msa_header_fname<<std::endl;
        exit(1);
    }
    msa_header_fp.close();
}

void MSA::load_graph_info(std::ifstream& stream) {
    std::ios::sync_with_stdio(false);
    std::string line;

    // first get length
    std::getline(stream,line);
    this->msa_len = std::stoi(line);

    // next get number of references
    std::getline(stream,line);
    this->num_refs = std::stoi(line);

    // next get the command used for building the graph
    std::getline(stream,line);
    this->build_cmd = line;

    this->graph = MSA_Graph(this->msa_len,this->num_refs); // initialize the graph here
}

void MSA::load_graph_contig_info(std::ifstream& stream) {
    // after the graph has been initialized we can now proceed to provide contig information
    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss("");
    std::string ref_id, ref_name,positions;
    while(std::getline(stream,line)) { // iterate over references
        ss.str(line);
        ss.clear();

        std::getline(ss, ref_name, '\t');
        std::getline(ss, ref_id, '\t');
        std::getline(ss,positions,'\t');
        this->graph.add_ref(ref_name,std::stoi(ref_id));

        std::stringstream pos_ss(positions);
        std::string pos;
        while (std::getline(pos_ss, pos, ',')) {
            this->graph.add_pos(std::stoi(ref_id),std::stoi(pos));
        }
    }
}

void MSA::_load_graph(std::ifstream& stream) {
    // now that everything else has been initialized - time to load vertices and edges into the graph
    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss(""),sub_ss("");
    std::string vt_pos, vt_contents, ref_id, nt, edge;
    MSA_Vertex mv;

    while(std::getline(stream,line)) { // iterate over references
        ss.str(line);
        ss.clear();

        std::getline(ss, vt_pos, '\t');
        std::getline(ss, vt_contents, '\n');

        ss.str(vt_contents);
        ss.clear();

        mv = MSA_Vertex(this->num_refs,std::stoi(vt_pos));

        while(std::getline(ss,vt_contents,'\t')) { // iterate over all contents of the vertex
            sub_ss.str(vt_contents);
            sub_ss.clear();
            std::getline(sub_ss, ref_id, ':');
            std::getline(sub_ss, nt, ';');
            std::getline(sub_ss, edge);

            mv.add_snp(nt,std::stoi(ref_id));
            if(std::stoi(edge) != 0){
                mv.add_edge(std::stoi(edge),std::stoi(ref_id));
            }
        }
        this->graph.add_vertex(std::stoi(vt_pos),mv);
    }

    this->graph.init_refcouts();
}

bool MSA::isMod(bam1_t* in_rec){
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS || opcode==BAM_CDEL || opcode==BAM_CREF_SKIP){
            return true;
        }
        else{
            continue;
        }
    }
    return false;
}

void MSA::write_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header){
    std::string ref_name = std::string(in_al_hdr->target_name[in_rec->core.tid]);
    in_rec->core.tid = 0;
    in_rec->core.pos = this->graph.get_new_position(ref_name,in_rec->core.pos);

    int ret_val = sam_write1(outSAM, outSAM_header, in_rec);
}

bool MSA::change_data(bam1_t *in_rec,int num_cigars,int* cigars,int cur_start,int cur_len,bool shift){
    int old_num_cigars = in_rec->core.n_cigar;

    int first_byte = cur_start/2;
    shift=false;
    if(cur_start%2==1){
        shift=true;
    }
    int num_bytes = (cur_len+1)/2;
    if(shift){
        num_bytes++;
    }

    int res_num_bytes = cur_len/2;
    if(cur_len%2==1){
        res_num_bytes++;
    }

    int data_len = (in_rec)->core.l_qname + (num_cigars * 4) + res_num_bytes + cur_len + (in_rec->l_data - (((old_num_cigars*4) + (in_rec)->core.l_qname + (((in_rec)->core.l_qseq + 1)/2) + (in_rec)->core.l_qseq)));

    int m_data=std::max(data_len,(int)in_rec->m_data);
    kroundup32(m_data);

    auto* data = (uint8_t*)calloc(m_data,1);

    // copy everything (QNAME) until CIGAR data
    int copy1_len = (uint8_t*)bam_get_cigar(in_rec) - in_rec->data;
    memcpy(data, in_rec->data, copy1_len);

    // copy CIGAR data
    int copy2_len = num_cigars * 4;
    memcpy(data + copy1_len, cigars, copy2_len);

    // shorten sequence string and copy it over to the new data
    uint8_t seq_ar[num_bytes];
    memcpy(seq_ar,bam_get_seq(in_rec)+first_byte,num_bytes);
    uint8_t *seq = seq_ar;

    if(shift){
        while (seq < seq_ar+num_bytes-1) {
            *seq = (*(seq)&0x0F)<<4 | (*(seq+1)&0xF0)>>4; // 0000 1111 - 1111 0000
            seq++;
        }
        *seq = (*(seq)&0x0F)<<4; // 0000 1111
    }

    memcpy(data + copy1_len + copy2_len, seq_ar, res_num_bytes);

    // shorten quality string and copy it over to the new data
    memcpy(data + copy1_len + copy2_len + res_num_bytes, bam_get_qual(in_rec)+cur_start, cur_len);

    // copy the remainder of the source data to the destination
    int copied_len = copy1_len + (old_num_cigars * 4) + ((in_rec->core.l_qseq + 1)>>1) + in_rec->core.l_qseq;
    int remain_len = in_rec->l_data - copied_len;
    memcpy(data + copy1_len + copy2_len + res_num_bytes + cur_len, in_rec->data+copied_len, remain_len);

    in_rec->core.n_cigar = num_cigars;

    free(in_rec->data);
    in_rec->data = data;
    in_rec->l_data = data_len;
    in_rec->m_data = m_data;
    in_rec->core.l_qseq = cur_len;

    return shift^((cur_len%2)==1);
}

void MSA::add_orig_ref_tags(bam1_t* in_rec,int ref,int new_end){
    uint8_t* ptr_op=bam_aux_get(in_rec,"ZA");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZA",'i',4,(uint8_t*)&ref);

    ptr_op=bam_aux_get(in_rec,"ZB");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZB",'i',4,(uint8_t*)&new_end);
}

void MSA::add_split_tags(bam1_t* in_rec,int cur_slice,int opcode,int length){
    uint8_t* ptr_op=bam_aux_get(in_rec,"ZC");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZC",'i',4,(uint8_t*)&cur_slice);

    ptr_op=bam_aux_get(in_rec,"ZD");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZD",'i',4,(uint8_t*)&opcode);

    ptr_op=bam_aux_get(in_rec,"ZJ");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZJ",'i',4,(uint8_t*)&length);
}

void MSA::add_seq_slice_tag(bam1_t* in_rec,uint8_t* seq_slice,int seq_len){
    uint8_t* ptr_op=bam_aux_get(in_rec,"ZL");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZL",'i',4,(uint8_t*)&seq_len);

    ptr_op=bam_aux_get(in_rec,"ZS");
    if(ptr_op){bam_aux_del(in_rec,ptr_op);}
    bam_aux_append(in_rec,"ZS",'Z',seq_len,seq_slice);
}

bool MSA::get_seq_slice(bam1_t* in_rec,bam1_t* out_rec,int cur_start,int cur_len,bool shift){
    int first_byte = cur_start/2;
    shift=false;
    if(cur_start%2==1){
        shift=true;
    }
    int num_bytes = (cur_len+1)/2;
    if(shift){
        num_bytes++;
    }

    int res_num_bytes = cur_len/2;
    if(cur_len%2==1){
        res_num_bytes++;
    }

    // shorten sequence string and copy it over to the new data
    uint8_t seq_ar[num_bytes];
    memcpy(seq_ar,bam_get_seq(in_rec)+first_byte,num_bytes);
    uint8_t *seq = seq_ar;

    if(shift){
        while (seq < seq_ar+num_bytes-1) {
            *seq = (*(seq)&0x0F)<<4 | (*(seq+1)&0xF0)>>4; // 0000 1111 - 1111 0000
            seq++;
        }
        *seq = (*(seq)&0x0F)<<4; // 0000 1111
    }
    seq=seq_ar;
    seq[res_num_bytes] = '\0';

    uint8_t* ptr_op=bam_aux_get(out_rec,"ZS");
    if(ptr_op){bam_aux_del(out_rec,ptr_op);}
    bam_aux_append(out_rec,"ZS",'Z',res_num_bytes+1,seq);

    // now also get the quality string into another aux tag
    uint8_t qual_ar[cur_len+1];
    memcpy(qual_ar,bam_get_qual(in_rec)+cur_start,cur_len);
    uint8_t* qual=qual_ar;
    qual[cur_len]='\0';
    ptr_op=bam_aux_get(out_rec,"ZQ");
    if(ptr_op){bam_aux_del(out_rec,ptr_op);}
    bam_aux_append(out_rec,"ZQ",'Z',cur_len+1,qual);

    return shift^((cur_len%2)==1);
}

void MSA::remove_aux_tags(bam1_t* rec){
    uint8_t* ptr_op=bam_aux_get(rec,"ZS");
    if(ptr_op){bam_aux_del(rec,ptr_op);}

    ptr_op=bam_aux_get(rec,"ZC");
    if(ptr_op){bam_aux_del(rec,ptr_op);}

    ptr_op=bam_aux_get(rec,"ZD");
    if(ptr_op){bam_aux_del(rec,ptr_op);}

    ptr_op=bam_aux_get(rec,"ZJ");
    if(ptr_op){bam_aux_del(rec,ptr_op);}
}

// split read based on the cigar string (Deletion/Insertion/SpliceSite)
void MSA::split_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header){
    int cigars[MAX_CIGARS];
    int cur_total_pos = in_rec->core.pos; // same as cur_pos but includes the soft clipping bases
    int cur_start=in_rec->core.pos,cur_end=in_rec->core.pos;
    int cur_local_start = 0;
    int cur_len = 0;
    int cur_ref_len=0;
    int num_cigars = 0;

    bam1_t* new_rec = bam_init1();
    new_rec = bam_dup1(in_rec);

    bool shift = false;
    int cur_slice = 0; // which slice within read is this
    int opcode;

    int new_end,new_start;

    std::string ref_name = std::string(in_al_hdr->target_name[in_rec->core.tid]);
    int refID = this->graph.get_id(ref_name); // reference ID of the read

    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS || opcode==BAM_CDEL || opcode==BAM_CREF_SKIP){
            shift = change_data(new_rec,num_cigars,cigars,cur_local_start,cur_len,shift);
            new_rec->core.pos = cur_start;
            new_rec->core.l_qseq = cur_len;
            new_end = this->graph.get_new_position(refID,bam_endpos(new_rec)-1); // bam_end_pos returns the coordinate of first base after the alignment which is why we subtract 1 from the returned value
            new_start = this->graph.get_new_position(refID,new_rec->core.pos);
            add_orig_ref_tags(new_rec,refID,new_end);
            add_split_tags(new_rec,cur_slice,opcode,length);
            // need to know where we are on the read and then use that in order to get the sequence of the insertion and correctly position the sliced fragment on the reference graph
            cur_local_start += cur_len;
            cur_start=cur_end;
            if(opcode==BAM_CINS){ // add the slice of the sequence here which can later be reinserted back into the joined read
                shift = get_seq_slice(in_rec,new_rec,cur_local_start,length,shift);
                cur_local_start+=length;
            }
            else{
                cur_end+=length;
                cur_start+=length;
            }
            new_rec->core.pos = new_start;
            new_rec->core.tid = 0;
            int ret_val = sam_write1(outSAM, outSAM_header, new_rec);
            cur_slice++;

            new_rec = bam_dup1(in_rec);
            num_cigars = 0;
            cur_len = 0;
            cur_ref_len=0;
            if(opcode!=BAM_CINS){
                cur_ref_len+=length;
            }
        }
        else{
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            cur_len+=length;
            if(opcode!=BAM_CSOFT_CLIP){
                cur_ref_len+=length;
                cur_end+=length;
            }
        }
    }
    shift = change_data(new_rec,num_cigars,cigars,cur_local_start,cur_len,shift);
    new_rec->core.pos = cur_start;
    new_rec->core.l_qseq = cur_len;
    new_end = this->graph.get_new_position(refID,bam_endpos(new_rec)-1); // bam_end_pos returns the coordinate of first base after the alignment which is why we subtract 1 from the returned value
    new_start = this->graph.get_new_position(refID,new_rec->core.pos);
    add_orig_ref_tags(new_rec,refID,new_end);
    new_rec->core.pos = new_start;
    new_rec->core.tid = 0;
    add_split_tags(new_rec,cur_slice,0,0); // 0 opcode here means that this is the last segment in the read
    int ret_val = sam_write1(outSAM, outSAM_header, new_rec);
}

// clean the graph based on the specified position
// TODO: not significant at the moment multiple strategies can be employed:
//    1. remove everything before current start
//    2. remove everything before current position but keep preceeding bases from the same reference as the first read
//    3. remove everything before current position but keep preceeding bases from the specified reference
//    4. do not remove anything
void MSA::clean(){
    // find first mapped position
    int first_pos_read=0,first_refID;
    this->graph.get_first_mapped_pos(first_pos_read,first_refID);

    int first_pos = 0;
    // 1. find the reference ID of the first mapping
    int first_pos_default = 0;
    if(this->gap_fillID==-1){
        first_pos_default = this->graph.get_first_pos(first_refID);
        first_pos = first_pos_default;
    }

    // 2. find where that reference begins in the graph
    //     a. if gapfillname is set
    int first_pos_ref = 0;
    if(this->gap_fillID!=-1){
        first_pos_ref = this->graph.get_first_pos(this->gap_fillID);
        if(first_pos_read<first_pos_ref){ // read mapping preceeds the requested genome - the assembly will begin directly with the read
            first_pos = first_pos_read;
        }
        else{
            first_pos = first_pos_ref;
        }
    }

    // 3. remove anything that preceeds the identified position
    this->graph.set_removed(0,first_pos);

    // now to clean the end sequence
    int last_pos;
    int last_pos_read,last_refID_read;
    this->graph.get_last_mapped_pos(last_pos_read,last_refID_read);

    int last_pos_ref = 0;
    bool set_by_gapfill = false;
    if(this->gap_fillID!=-1){
        set_by_gapfill = true;
        first_pos_ref = this->graph.get_last_pos(this->gap_fillID);
        if(last_pos_read<last_pos_ref){ // read mapping preceeds the requested genome - the assembly will begin directly with the read
            last_pos = last_pos_read;
        }
        else{
            last_pos = last_pos_ref;
        }
    }
    if(!set_by_gapfill){
        last_pos = last_pos_read;
    }
    this->graph.set_removed(last_pos,this->graph.get_len()-1);

    // 4. for each gap and ends do the same process
}

void MSA::change_cigar(bam1_t* in_rec,int s){

}

void MSA::l2range(std::vector<int>& l,std::vector<std::pair<int,int>>& r){
    std::sort(l.begin(),l.end());
    r.emplace_back(std::make_pair(l.front(),l.front()));
    bool new_entry = false;
    for(int i=0;i<l.size();i++){
        if(new_entry){
            r.emplace_back(std::make_pair(l[i],l[i]));
            new_entry = false;
            continue;
        }
        if(i==l.size()){
            return;
        }
        else{
            if(l[i+1]-l[i]==1){
                r.back().second++;
            }
            else if(l[i+1]-l[i]<=0){
                continue;
            }
            else{
                new_entry=true;
            }
        }
    }
}

void MSA::add_cigar(bam1_t *curAl,int num_cigars,int* cigars){
    int old_num_cigars = curAl->core.n_cigar;
    int data_len=curAl->l_data+4*(num_cigars-old_num_cigars);
    int m_data=std::max(data_len,(int)curAl->m_data);
    kroundup32(m_data);

    auto* data = (uint8_t*)calloc(m_data,1);

    int copy1_len = (uint8_t*)bam_get_cigar(curAl) - curAl->data;
    memcpy(data, curAl->data, copy1_len);

    int copy2_len = num_cigars * 4;
    memcpy(data + copy1_len, cigars, copy2_len);

    int copy3_len = curAl->l_data - copy1_len - (old_num_cigars * 4);
    memcpy(data + copy1_len + copy2_len, bam_get_seq(curAl), copy3_len);

    curAl->core.n_cigar = num_cigars;

    free(curAl->data);
    curAl->data = data;
    curAl->l_data = data_len;
    curAl->m_data = m_data;
}

// TODO: when filling in the gaps in the consensus what we can do by default is:
//    check which genome has the highest number of unique mappings and use that

void MSA::create_del(bam1_t* in_rec,std::vector<std::pair<int,int>>& not_removed){
    // currently we assume that positions in not_removed are sorted (which should hold true, since the alignment is evaluated in order)
    // however, this may lead to errors

    int cigars[MAX_CIGARS];
    int num_cigars=0;

    int nr_idx = 0;

//    int cur_ref_pos=in_rec->core.pos; // same as cur_pos but includes the soft clipping bases
    int cur_read_pos=0;
    int cur_added_pos=0; // tracks the numbe rof bases which have already been moved from the old cigar into the new one
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if(opcode==BAM_CSOFT_CLIP || opcode==BAM_CMATCH || opcode==BAM_CINS){
                cur_read_pos+=length;
        }
        int pre_len=0,post_len=0;
        bool was_set=false;
        while(cur_read_pos>not_removed[nr_idx].first && nr_idx!=not_removed.size()){ // add all the not_removed bases as relevant deletions
            was_set=true;
            pre_len = length - (cur_read_pos-not_removed[nr_idx].first)-cur_added_pos;
            post_len = length - pre_len - cur_added_pos;
            cigars[num_cigars]=opcode|(pre_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            if(opcode==BAM_CSOFT_CLIP || opcode==BAM_CMATCH || opcode==BAM_CINS){
                // consumes both reference and query
                cur_added_pos+=pre_len;
            }

            int nr_len = (not_removed[nr_idx].second-not_removed[nr_idx].first)+1;
            cigars[num_cigars]=BAM_CDEL|(nr_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;

//            cur_read_pos+=nr_len;

            nr_idx++;
        }
        if(was_set){
            cigars[num_cigars]=opcode|(post_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            was_set=false;
        }
        else{
            if(cur_read_pos<=not_removed[nr_idx].first){
                cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
                ++num_cigars;
            }
        }
    }

    add_cigar(in_rec,num_cigars,cigars);
    return;
}

void MSA::create_ins_old(bam1_t* in_rec,std::vector<std::pair<int,int>>& added){
    int cigars[MAX_CIGARS];
    int num_cigars=0;

    int nr_idx = 0;
    int cur_read_pos=0; // delayed value does never advances ahead
    int last_read_pos = 0; // the total number of bases in the cigar operations fully consumed thus far
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
//        if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){std::cout<<"p1"<<std::endl;}
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if(opcode==BAM_CSOFT_CLIP || opcode==BAM_CMATCH || opcode==BAM_CINS){
            // consumes both reference and query
            cur_read_pos+=length;
        }
        int pre_len=0,post_len=length,nr_len=0,cum_pre_len=0;
        bool was_set=false;
        bool ins_created=false;
        bool skip_post=false;
        while(cur_read_pos>=added[nr_idx].first && nr_idx!=added.size()){ // add all the not_removed bases as relevant deletions
//            if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){std::cout<<"p2"<<std::endl;}
            was_set=true;
            skip_post=false;
            pre_len = length - ((cur_read_pos-added[nr_idx].first)+cum_pre_len);
            cum_pre_len+=pre_len;
            if(ins_created){
                post_len-=nr_len;
                ins_created=false;
            }
            cigars[num_cigars]=opcode|(pre_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;

            nr_len = (added[nr_idx].second-added[nr_idx].first)+1;
            cum_pre_len+=nr_len;
            // since insertions advance position on the read we need to make sure the modification will not step over the bounds of the read and adjust accordingly
            if(pre_len+nr_len>length){
                nr_len = nr_len - ((pre_len+nr_len)-length);
                skip_post = true;
            }
            cigars[num_cigars]=BAM_CINS|(nr_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            ins_created=true;
            nr_idx++;
        }
        if(was_set && !skip_post){
//            if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){std::cout<<"p3"<<std::endl;}
            if(ins_created){
                post_len-=nr_len;
                ins_created=false;
            }
            post_len = length - cum_pre_len;
            cigars[num_cigars]=opcode|(post_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            was_set=false;
        }
        else if(skip_post){
//            if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){std::cout<<"p4"<<std::endl;}
            break;
        }
        else{
            if(cur_read_pos<added[nr_idx].first){
//                if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){std::cout<<"p5"<<std::endl;}
                cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
                ++num_cigars;
            }
        }
    }

    add_cigar(in_rec,num_cigars,cigars);
    return;
}

void MSA::create_ins(bam1_t* in_rec,std::unordered_set<int>& added){
    int cigars[MAX_CIGARS];
    int num_cigars=0;

    std::vector<uint8_t> blown_up_cigars;

    // first need to populate the blownup containers
    uint32_t *cigar_full=bam_get_cigar(in_rec);
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
        std::fill_n (std::back_inserter(blown_up_cigars), length, opcode);
    }

    int prev_op = -1; // previous operation
    int ins_pos_inc=0; // for deletions and introns tells by how much the position needs to advance for insertions

    int cur_len = 0; // length of the current operation

    int op;
    for(int i=0;i<blown_up_cigars.size();i++){
        op = blown_up_cigars[i];
        if(added.find(i-ins_pos_inc)!=added.end()){
            if(op == BAM_CDEL || op == BAM_CREF_SKIP){
                ins_pos_inc++;
            }
            else{
                op = BAM_CINS;
            }
        }
        if(op!=prev_op){
            cur_len=1;
            cigars[num_cigars]=op|(cur_len<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            prev_op=op;
        }
        else{
            cur_len++;
            cigars[num_cigars-1]=op|(cur_len<<BAM_CIGAR_SHIFT);
        }
    }

    add_cigar(in_rec,num_cigars,cigars);
    return;
}

// parse fitted read and adjust the graph
void MSA::parse_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header){
    uint8_t* ptr_nm_1=bam_aux_get(in_rec,"ZB");
    int tag_ref_end = bam_aux2i(ptr_nm_1);
    ptr_nm_1=bam_aux_get(in_rec,"ZA");
    int tag_refID = bam_aux2i(ptr_nm_1);

    int in_rec_ref_start = in_rec->core.pos;

    int new_start,s;
    std::vector<int> not_removed_tmp;
    std::unordered_set<int> added_tmp;
    std::vector<std::pair<int,int>> not_removed,added;

//    if(std::strcmp(bam_get_qname(in_rec),"AF084936_9472_9622_0:1:0_0:1:0_e8a/1")==0){
//        std::cout<<"found"<<std::endl;
//    }

    this->graph.fit_read(tag_refID,in_rec_ref_start,tag_ref_end,new_start,s,not_removed_tmp,added_tmp);
    in_rec_ref_start = new_start;

    if(in_rec_ref_start != NULL){
        in_rec->core.pos = in_rec_ref_start;
        if(s>0){
            change_cigar(in_rec,s);
        }
        if(not_removed_tmp.size()>0){
            // instead of introducing a deletion into the actual cigar string - perhaps would make sense to simply split the read and record the type for later
            // how can we do this without having to repeat the same process for the reads that fall in here
            l2range(not_removed_tmp,not_removed);
//            if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){
//                print_cigar(in_rec);
//            }
            create_del(in_rec,not_removed);
//            if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){
//                print_cigar(in_rec);
//            }
        }
        if(added_tmp.size()>0){
            // instead of introducing a deletion into the actual cigar string - perhaps would make sense to simply split the read and record the type for later
//            l2range(added_tmp,added);
            create_ins(in_rec,added_tmp);
        }
//        if(std::strcmp(bam_get_qname(in_rec),"KF234628_7421_7571_0:0:0_0:0:0_78c/1")==0){
//            print_cigar(in_rec);
//        }
        int ret_val = sam_write1(outSAM, outSAM_header, in_rec);
    }
}

void MSA::join_cigars_old(std::vector<bam1_t*>& reads,uint8_t *cigars,int& new_n_cigar_bytes){
    int cur_mem_pos = 0;

    for(int i=0;i<reads.size();i++){
        bam1_t* rec = reads[i];
        int cigar_len = rec->core.n_cigar * 4;
        uint32_t *cigar_full=bam_get_cigar(rec);
        int opcode=bam_cigar_op(cigar_full[0]);
        int length=bam_cigar_oplen(cigar_full[0]);

        if(i==reads.size()-1){ // last fragment
            if(opcode==BAM_CDEL){ // deletions will be computed based on the distance to the nearest start position
                memcpy(cigars + cur_mem_pos, &cigar_full[1], cigar_len-4);
                cur_mem_pos+=(cigar_len-4);
                new_n_cigar_bytes+=(cigar_len-4);
            }
            else{
                memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
                cur_mem_pos+=cigar_len;
                new_n_cigar_bytes+=cigar_len;
            }
            break;
        }
        if(opcode==BAM_CINS){
            bam1_t* next_read = reads[i+1];
            int next_start=next_read->core.pos;
            int cur_end=bam_endpos(rec);
            int nd_len=next_start-(cur_end+1);
            memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
            cur_mem_pos+=cigar_len;
            new_n_cigar_bytes+=cigar_len;

            if(nd_len>0){
                cigar_full[0] = BAM_CDEL|(nd_len<<BAM_CIGAR_SHIFT);
                memcpy(cigars + cur_mem_pos, cigar_full, 4);
                cur_mem_pos+=4;
                new_n_cigar_bytes+=4;
            }
        }
        else{
            bam1_t* next_read = reads[i+1];
            int next_start=next_read->core.pos;
            int cur_end=bam_endpos(rec);
            int nd_len=next_start-cur_end;
            if(opcode==BAM_CDEL){ // deletions will be computed based on the distance to the nearest start position
                memcpy(cigars + cur_mem_pos, &cigar_full[1], cigar_len-4);
                cur_mem_pos+=(cigar_len-4);
                new_n_cigar_bytes+=(cigar_len-4);
            }
            else{
                memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
                cur_mem_pos+=cigar_len;
                new_n_cigar_bytes+=cigar_len;
            }
            if(nd_len>0) {
                cigar_full[0] = BAM_CDEL|(nd_len<<BAM_CIGAR_SHIFT);
                memcpy(cigars + cur_mem_pos, cigar_full, 4);
                cur_mem_pos+=4;
                new_n_cigar_bytes+=4;
            }
        }
    }
}

void MSA::get_split_op_tags(bam1_t* rec,int& opcode,int&oplen){
    uint8_t* ptr_op=bam_aux_get(rec,"ZD");
    if(ptr_op){
        opcode = bam_aux2i(ptr_op);
        if(opcode==1){ // insertion
            ptr_op=bam_aux_get(rec,"ZJ");
            if(ptr_op){
                oplen = bam_aux2i(ptr_op);
            }
            else{
                std::cerr<<"ERROR:: did not find corresponding length"<<std::endl;
                exit(-1);
            }
        }
    }
}

// compute what opcode and oplen are required to correctly join the two fragments
void MSA::get_dist2next(bam1_t* rec,bam1_t* next,int opcode,int oplen){
    oplen = next->core.pos - bam_endpos(rec);
    if(oplen == 0){ // no insertions or deletions needed - can simply write as 0M
        opcode = BAM_CMATCH;
        return;
    }
    else if(oplen>0){
        opcode = BAM_CDEL;
        return;
    }
    else{ // inferred insertion
        opcode = BAM_CINS;
        oplen = 0-oplen; // taking the absolute
        return;
    }
}

void MSA::join_cigars(std::vector<bam1_t*>& reads,uint8_t *cigars,int& new_n_cigar_bytes,std::unordered_set<int>& added){
    int cur_mem_pos = 0;

//    if(std::strcmp(bam_get_qname(reads[0]),"K03455_8848_8998_0:1:0_0:1:0_b0a/1")==0){
//        std::cout<<"found"<<std::endl;
//    }

    int read_pos=0; // position on the read

    for(int i=0;i<reads.size();i++){
        bam1_t* rec = reads[i];
        int cigar_len = rec->core.n_cigar * 4;
        uint32_t *cigar_full=bam_get_cigar(rec);
        int opcode,oplen;
        get_split_op_tags(rec,opcode,oplen);

        if(i==reads.size()-1){ // last fragment
            memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
            cur_mem_pos+=cigar_len;
            new_n_cigar_bytes+=cigar_len;
            break;
        }
        if(opcode==BAM_CINS){
            bam1_t* next_read = reads[i+1];
            int next_start=next_read->core.pos;
            int cur_end=bam_endpos(rec);
            int nd_len=next_start-cur_end;

            // copy the cigar of the current slice
            memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
            read_pos+=rec->core.l_qseq;
            cur_mem_pos+=cigar_len;
            new_n_cigar_bytes+=cigar_len;

            if(nd_len<0){ // need to extend the insertion due to the overlap
                // what if here we just remember positions to create insertions and deletions, and then introduce them using conventional algorithms used before this step when fitting the segments onto the graph?
                for(int cp=(read_pos+nd_len);cp<read_pos;cp++){
                    added.insert(cp);
                }
            }

            // create the insertion
            cigar_full[0] = BAM_CINS|(oplen<<BAM_CIGAR_SHIFT);
            read_pos+=oplen;
            memcpy(cigars + cur_mem_pos, cigar_full, 4);
            cur_mem_pos+=4;
            new_n_cigar_bytes+=4;

            if(nd_len>0){
                cigar_full[0] = BAM_CDEL|(nd_len<<BAM_CIGAR_SHIFT);
                memcpy(cigars + cur_mem_pos, cigar_full, 4);
                cur_mem_pos+=4;
                new_n_cigar_bytes+=4;
            }
        }
        else{
            bam1_t* next_read = reads[i+1];
            int next_start=next_read->core.pos;
            int cur_end=bam_endpos(rec);
            int nd_len=next_start-cur_end;

            // copy the cigar of the current slice
            memcpy(cigars + cur_mem_pos, cigar_full, cigar_len);
            read_pos+=rec->core.l_qseq;
            cur_mem_pos+=cigar_len;
            new_n_cigar_bytes+=cigar_len;

            if(nd_len>=0){
                cigar_full[0] = opcode|(nd_len<<BAM_CIGAR_SHIFT);
                memcpy(cigars + cur_mem_pos, cigar_full, 4);
                cur_mem_pos+=4;
                new_n_cigar_bytes+=4;
            }
        }
    }
}

int MSA::get_ins_len(bam1_t* rec){
    uint8_t* ptr_op=bam_aux_get(rec,"ZD");
    if(ptr_op){
        int opcode = bam_aux2i(ptr_op);
        if(opcode==1){ // insertion
            ptr_op=bam_aux_get(rec,"ZJ");
            if(ptr_op){
                int oplen = bam_aux2i(ptr_op);
                return oplen;
            }
        }
    }
    return 0;
}

void MSA::get_ins_seq(bam1_t* rec,uint8_t* seq,int& seq_len){
    uint8_t* ptr_seq=bam_aux_get(rec,"ZS"); // ZS is only set if insertion
    if(ptr_seq){
        uint8_t* ptr_oplen=bam_aux_get(rec,"ZJ");
        if(ptr_oplen){
            seq_len = bam_aux2i(ptr_oplen);
                memcpy(seq,bam_aux2Z(ptr_seq),(seq_len+1)/2);
        }
        else{
            std::cerr<<"ERROR:::sequence exists but no length"<<std::endl;
            exit(-1);
        }
    }
}

void MSA::merge_seqs(uint8_t* data,uint8_t* seq,int seq_len,int& cur_mem_pos,bool orphan){
    int num_bytes = (seq_len+1)/2;
    uint8_t* seq_dup = seq;
    if(orphan){
        uint8_t byte_ar[1];
        memcpy(byte_ar,data+cur_mem_pos-1,1); // copy the last byte over which will contain the orphan bit

        uint8_t *byte = byte_ar;
        *byte = (*(byte)&0xF0) | (*(seq)&0xF0)>>4;
        memcpy(data+cur_mem_pos-1,byte_ar,1);

        int i=0;
        while (i < num_bytes-1) {
            *seq = (*(seq)&0x0F)<<4 | (*(seq+1)&0xF0)>>4; // 0000 1111 - 1111 0000
            seq++;
            i++;
        }
        *seq = (*(seq)&0x0F)<<4; // 0000 1111
        seq=seq_dup;
    }

    int last_bit = int(orphan&(seq_len%2==1));
    memcpy(data+cur_mem_pos,seq,num_bytes-last_bit);
    cur_mem_pos+=(num_bytes-last_bit);
}

void MSA::join_seqs(std::vector<bam1_t*>& reads,uint8_t* data,int& cur_mem_pos,int max_seq_len){
    // iteratively add sequence pieces to the new record
    bool merge_last_first = false; // if set notifies that first 4 bytes of the last byte of previous sequence need to be combined accordingly with the first 4 bytes of the current sequence before copying the rest
    bam1_t *rec;
    for(int i=0;i<reads.size();i++){
        rec = bam_dup1(reads[i]);
        int seq_len = rec->core.l_qseq;

        int num_bytes = (seq_len+1)/2;
        uint8_t seq_ar[num_bytes];
        memcpy(seq_ar,bam_get_seq(rec),num_bytes);
        uint8_t *seq = seq_ar;

        if(merge_last_first){
            merge_seqs(data,seq,seq_len,cur_mem_pos,true);
        }
        else{
            memcpy(data+cur_mem_pos,seq_ar,num_bytes);
            cur_mem_pos+=num_bytes;
        }
        // process potential insertions here
        uint8_t seq_tag_ar[(max_seq_len+1)/2];
        uint8_t* seq_tag = seq_tag_ar;
        int seq_len_tag = 0;
        get_ins_seq(rec,seq_tag,seq_len_tag);
        if(seq_len_tag>0){ // insertion exists - need to add accordingly
            bool shift = merge_last_first;
            if(seq_len%2==1){
                shift=!shift;
            }
            merge_seqs(data,seq_tag,seq_len_tag,cur_mem_pos,shift);
            seq_len+=seq_len_tag; // append the length of the modification to the rest for proper evaluation of shifting below
        }
        if(seq_len%2==1){
            merge_last_first=!merge_last_first;
        }
    }
}

void MSA::get_ins_qual(bam1_t* rec,uint8_t* qual,int& qual_len){
    uint8_t* ptr_qual=bam_aux_get(rec,"ZQ"); // ZQ is only set if insertion
    if(ptr_qual){
        uint8_t* ptr_oplen=bam_aux_get(rec,"ZJ");
        if(ptr_oplen){
            qual_len = bam_aux2i(ptr_oplen);
            memcpy(qual,bam_aux2Z(ptr_qual),qual_len);
        }
        else{
            std::cerr<<"ERROR:::quality exists but no length"<<std::endl;
            exit(-1);
        }
    }
}

void MSA::join_quals(std::vector<bam1_t*>& reads,uint8_t* data,int& cur_mem_pos,int max_seq_len){
    bam1_t *rec;
    for(int i=0;i<reads.size();i++){
        rec = bam_dup1(reads[i]);
        int qual_len = rec->core.l_qseq;

        memcpy(data+cur_mem_pos,bam_get_qual(rec),qual_len);
        cur_mem_pos+=qual_len;
        uint8_t qual_tag_ar[max_seq_len];// = (uint8_t*)calloc(max_seq_len,1);
        uint8_t* qual_tag = qual_tag_ar;
        int qual_len_tag = 0;
        get_ins_qual(rec,qual_tag,qual_len_tag);
        if(qual_len_tag>0){ // insertion exists - need to add accordingly
            memcpy(data+cur_mem_pos,qual_tag,qual_len_tag);
            cur_mem_pos+=qual_len_tag;
        }
    }
}

void MSA::joinReads(std::vector<bam1_t*>& reads,samFile *outSAM_joined,bam_hdr_t *outSAM_joined_header){
    if(reads.empty()){
        return;
    }
    // get the size of each of the data blocks
    int total_seq_len = 0;
    int total_num_cigars = 0;
    for(bam1_t *rec : reads){
        int ins_len = get_ins_len(rec);
        total_seq_len+=rec->core.l_qseq+ins_len;
    }

    int seq_num_bytes = total_seq_len/2;
    if(total_seq_len%2==1){ // because there is a shift at both ends - the results will contain one byte less
        seq_num_bytes++;
    }

    // iteratively add cigar information here
    // need to treat Deletions and splice junctions as non-existent and replace their length with the actual length inside
    int cur_mem_pos = 0;
    int new_n_cigar_bytes = 0;

    uint8_t* cigars = (uint8_t*)calloc(MAX_CIGARS,1);

    std::unordered_set<int> added;
    join_cigars(reads,cigars,new_n_cigar_bytes,added);

    bam1_t* new_rec = bam_dup1(reads.front());

    int new_aux_len = new_rec->l_data - (((new_rec->core.n_cigar)*4) + new_rec->core.l_qname + ((new_rec->core.l_qseq + 1)/2) + new_rec->core.l_qseq);

    int data_len = new_rec->core.l_qname + new_n_cigar_bytes + seq_num_bytes + total_seq_len + new_aux_len; // the total length of the data to be populated from joined reads into the final record

    int m_data=std::max(data_len,(int)new_rec->m_data);
    kroundup32(m_data);

    auto* data = (uint8_t*)calloc(m_data,1);

    // copy everything (QNAME) until CIGAR data
    int name_len = (uint8_t*)bam_get_cigar(new_rec) - new_rec->data;
    memcpy(data, new_rec->data, name_len);
    cur_mem_pos = name_len;

    memcpy(data+cur_mem_pos, cigars, new_n_cigar_bytes);
    cur_mem_pos += new_n_cigar_bytes;

    join_seqs(reads,data,cur_mem_pos,total_seq_len);

    // TODO: need to take care of the quality strings and move them to the tags in case of insertions as well - should be fairly simple
    join_quals(reads,data,cur_mem_pos,total_seq_len);

    // lastly copy over the auxilary data
    int copied_len = name_len + (reads.front()->core.n_cigar * 4) + ((reads.front()->core.l_qseq+1)/2) + reads.front()->core.l_qseq;
    int remain_len = reads.front()->l_data - copied_len;
    memcpy(data+cur_mem_pos,reads.front()->data+copied_len,remain_len);

    new_rec->core.n_cigar = (new_n_cigar_bytes/4);
    new_rec->core.l_qseq = total_seq_len;

    free(new_rec->data);
    new_rec->data = data;
    new_rec->l_data = data_len;
    new_rec->m_data = m_data;

    remove_aux_tags(new_rec);

    if(!added.empty()){
        create_ins(new_rec,added);
    }

    int ret = sam_write1(outSAM_joined,outSAM_joined_header,new_rec);
    return;
}

void MSA::realign(std::string in_sam,std::string out_sam){
    samFile *msa_hdr_fp = hts_open(this->msa_header_fname.c_str(),"r");
    bam_hdr_t *msa_hdr = sam_hdr_read(msa_hdr_fp);

    std::cerr<<"@LOG::::Begin fitting alignments to MSA"<<std::endl;

    samFile *outSAM=sam_open(out_sam.c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_init();
    outSAM_header=bam_hdr_dup(msa_hdr);
    sam_hdr_write(outSAM,outSAM_header);

    samFile *in_al=sam_open(in_sam.c_str(),"r");
    bam_hdr_t *in_al_hdr = sam_hdr_read(in_al); //read header
    in_al_hdr->ignore_sam_err=1;
    bam1_t *in_rec = bam_init1(); //initialize an alignment
    int ret;

    while(sam_read1(in_al, in_al_hdr, in_rec) >= 0) {
        if(isMod(in_rec)){
            split_read(in_rec,in_al_hdr,outSAM,outSAM_header);
        }
        else{
            std::string ref_name = std::string(in_al_hdr->target_name[in_rec->core.tid]);
            int refID = this->graph.get_id(ref_name); // reference ID of the read
            int new_end = this->graph.get_new_position(refID,bam_endpos(in_rec)-1);  // bam_end_pos returns the coordinate of first base after the alignment which is why we subtract 1 from the returned value
            int new_start = this->graph.get_new_position(refID,in_rec->core.pos);
            add_orig_ref_tags(in_rec,refID,new_end);
            in_rec->core.pos = new_start;
            in_rec->core.tid = 0;
            int ret_val = sam_write1(outSAM, outSAM_header, in_rec);
        }
    }

    bam_destroy1(in_rec);
    sam_close(in_al);
    sam_close(outSAM);
    bam_hdr_destroy(outSAM_header);

    std::cerr<<"@LOG::::Done fitting alignment to MSA"<<std::endl;

    std::cerr<<"@LOG::::Begin sorting fitted alignment"<<std::endl;
    std::string sam_sort_cmd = "samtools sort -o "+out_sam+".sorted.bam "+out_sam;
    int res_sam_sort = system(sam_sort_cmd.c_str());
    std::cerr<<"@LOG::::Done sorting fitted alignment"<<std::endl;

    std::cerr<<"@LOG::::Begin cleaning graph"<<std::endl;
    std::string out_sam_clean = out_sam+".clean";
    samFile *outSAM_clean=sam_open(out_sam_clean.c_str(),"wb");
    bam_hdr_t *outSAM_clean_header=bam_hdr_init();
    outSAM_clean_header=bam_hdr_dup(msa_hdr);
    sam_hdr_write(outSAM_clean,outSAM_clean_header);
    bam_hdr_destroy(outSAM_clean_header);

    std::string sorted_al_fname = out_sam+".sorted.bam";
    in_al=sam_open(sorted_al_fname.c_str(),"r");
    in_al_hdr = sam_hdr_read(in_al); //read header
    in_al_hdr->ignore_sam_err=1;
    in_rec = bam_init1(); //initialize an alignment

//    this->clean(); // TODO: clean after parsing and then perform last correction after when joining
//                         this can be achieved by creating another vector to keep track of the vertices removed during the cleanup
//                         the modifications should not disturb the actual alignments

    while(sam_read1(in_al, in_al_hdr, in_rec) >= 0) {
//        std::cout<<bam_get_qname(in_rec)<<std::endl;
        parse_read(in_rec,in_al_hdr,outSAM_clean,outSAM_clean_header);
    }

    std::string consensus_fa_fname(out_sam);
    consensus_fa_fname.append(".cons.fasta");
    this->graph.save_merged_fasta(consensus_fa_fname);

    bam_destroy1(in_rec);
    sam_close(in_al);
    sam_close(outSAM_clean);
    std::cerr<<"@LOG::::Done cleaning graph"<<std::endl;

    std::cerr<<"@LOG::::Begin sorting disjoint alignments by name"<<std::endl;

    sam_sort_cmd = "samtools sort -n -o "+out_sam+".sorted_name.bam "+out_sam+".clean";
    res_sam_sort = system(sam_sort_cmd.c_str());

    std::cerr<<"@LOG::::Done sorting disjoint alignments by name"<<std::endl;

    std::cerr<<"@LOG::::Begin joining reads"<<std::endl;
    std::string out_sam_joined = out_sam+".joined";
    samFile *outSAM_joined=sam_open(out_sam_joined.c_str(),"wb");
    bam_hdr_t *outSAM_joined_header=bam_hdr_init();
    outSAM_joined_header=bam_hdr_dup(msa_hdr);
    sam_hdr_write(outSAM_joined,outSAM_joined_header);
    bam_hdr_destroy(outSAM_joined_header);

    sorted_al_fname = out_sam+".sorted_name.bam";
    in_al=sam_open(sorted_al_fname.c_str(),"r");
    in_al_hdr = sam_hdr_read(in_al); //read header
    in_al_hdr->ignore_sam_err=1;
    in_rec = bam_init1(); //initialize an alignment

    // join reads here and hope for the best
    std::string last_read;
    std::vector<bam1_t*> reads;

    while(sam_read1(in_al, in_al_hdr, in_rec) >= 0) {
        uint8_t* ptr_zc_1=bam_aux_get(in_rec,"ZC");
        if(!ptr_zc_1){ // the read was not split
            int ret_val = sam_write1(outSAM_joined, outSAM_joined_header, in_rec);
            continue;
        }
        else{
            // collect chains of fragments and join them according to the graph structure
            if(std::strcmp(bam_get_qname(in_rec),last_read.c_str())==0){
                reads.emplace_back(bam_dup1(in_rec));
            }
            else{ // can join reads together
                joinReads(reads,outSAM_joined,outSAM_joined_header);
                last_read = bam_get_qname(in_rec);
                reads.clear();
                reads.emplace_back(bam_dup1(in_rec));
            }
        }
    }

    bam_destroy1(in_rec);
    sam_close(in_al);
    sam_close(outSAM_joined);

    std::cerr<<"@LOG::::Done joining reads"<<std::endl;
}

void MSA::fit_annotation(std::string in_gff, std::string out_gff){
    std::cerr<<"@LOG::::Begin fitting annotation"<<std::endl;
    this->graph.fit_annotation(in_gff,out_gff);
    std::cerr<<"@LOG::""Done fitting annotation"<<std::endl;
    return;
}

void MSA::fit_bed(std::string in_bed, std::string out_bed){ // TODO: need BED parser
//    std::cerr<<"@LOG::::Begin fitting annotation"<<std::endl;
//    this->graph.fit_annotation(in_bed,out_bed);
//    std::cerr<<"@LOG::""Done loading annotation"<<std::endl;
    return;
}

void MSA::set_gapfillname(std::string ref_name){
    this->gap_fillID=this->graph.get_id(ref_name);
}

