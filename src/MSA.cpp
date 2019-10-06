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
        std::cerr<<"FATAL: Couldn't open graph contents file: "<<msa_header_fname<<std::endl;
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
        int old_pos = 0;
        while (std::getline(pos_ss, pos, ',')) {
            this->graph.add_pos(std::stoi(ref_id),old_pos,std::stoi(pos));
            old_pos++;
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

void print_cigar(bam1_t *al){
    for (uint8_t c=0;c<al->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(al);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
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
//    int32_t qlen = new_rec->core.l_qseq;
//    int8_t *buf = NULL;
//    buf = static_cast<int8_t *>(realloc(buf, qlen+1));
//    buf[qlen] = '\0';
//    uint8_t* seq = bam_get_qual(new_rec);
//    for (int i = 0; i < qlen; ++i)
//        buf[i] = bam_seqi(seq, i);
//    for (int i = 0; i < qlen; ++i) {
//        buf[i] = seq_nt16_str[buf[i]];
//    }
//    std::string str_seq((char*)(char*)buf);
    std::cout<<bam_get_qual(new_rec)<<std::endl;
}

bool MSA::change_data(bam1_t *in_rec,int num_cigars,int* cigars,int cur_start,int cur_len,bool shift){
    int old_num_cigars = in_rec->core.n_cigar;

    int data_len=in_rec->l_data;
    data_len = data_len - 4*(old_num_cigars - num_cigars); // modify with respect to the shortened new cigar string length
    data_len = data_len - (in_rec->core.l_qseq - cur_len)/2; // modify with respect to the shortened new sequence string
    data_len = data_len - (in_rec->core.l_qseq - cur_len); // modify with respect to the shortened new quality string

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
    int first_byte = cur_start/2;
    shift=false;
    if(cur_start%2==1){
        shift=true;
    }
    int num_bytes = cur_len/2;
    if(shift){
        num_bytes++;
    }
    if(cur_len%2==0){
        num_bytes++;
    }

    int res_num_bytes = cur_len/2;
    if(cur_len%2==1){ // becase there is a shift at both ends - the results will contain one byte less
        res_num_bytes++;
    }

//    if(std::strcmp(bam_get_qname(in_rec),"AF004885_119_269_0:0:2_0:0:2_78/1")==0) {
//        std::cout<<cur_start<<"\t"<<cur_len<<"\t"<<first_byte<<"\t"<<num_bytes<<"\t"<<res_num_bytes<<std::endl;
//    }

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
    int copied_len = copy1_len + (old_num_cigars * 4) + ((in_rec->core.l_qseq+1)/2) + in_rec->core.l_qseq;
    int remain_len = in_rec->l_data - copied_len;
    memcpy(data + copy1_len + copy2_len + res_num_bytes + cur_len, in_rec->data+copied_len, remain_len);

    in_rec->core.n_cigar = num_cigars;

    free(in_rec->data);
    in_rec->data = data;
    in_rec->l_data = data_len;
    in_rec->m_data = m_data;

    return shift^((cur_len%2)==1);
}

// split read based on the cigar string (Deletion/Insertion/SpliceSite)
void MSA::split_read(bam1_t* in_rec,bam_hdr_t *in_al_hdr,samFile* outSAM,bam_hdr_t* outSAM_header){
    int cigars[MAX_CIGARS];
    int cur_total_pos = in_rec->core.pos; // same as cur_pos but includes the soft clipping bases
    int cur_start=in_rec->core.pos;
    int cur_local_start = 0;
    int cur_len = 0;
    int num_cigars = 0;

    bam1_t* new_rec = bam_init1();
    new_rec = bam_dup1(in_rec);

    bool shift = false;
    for (uint8_t c=0;c<in_rec->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(in_rec);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS || opcode==BAM_CDEL || opcode==BAM_CREF_SKIP){
            shift = change_data(new_rec,num_cigars,cigars,cur_local_start,cur_len,shift);
            new_rec->core.pos = cur_start;
            new_rec->core.l_qseq = cur_len;
            write_read(new_rec,in_al_hdr,outSAM,outSAM_header);

//            if(std::strcmp(bam_get_qname(in_rec),"AF004885_119_269_0:0:2_0:0:2_78/1")==0) {
//                std::cout<<bam_get_qname(in_rec)<<std::endl;
//                std::cout<<cur_len<<std::endl;
//                std::cout<<cur_local_start<<std::endl;
//
//                print_seq(in_rec);
//                print_seq(new_rec);
//                print_cigar(in_rec);
//                print_cigar(new_rec);
//            }

            new_rec = bam_dup1(in_rec);
            num_cigars = 0;
            cur_start+=cur_len;
            cur_local_start += cur_len;
            cur_len = 0;
            if(opcode==BAM_CINS){
                cur_len+=length;
            }
            else{
                cur_start+=length;
            }
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
        }
        else{
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            cur_len+=length;
        }
    }
    shift = change_data(new_rec,num_cigars,cigars,cur_local_start,cur_len,shift);
    new_rec->core.pos = cur_start;
    new_rec->core.l_qseq = cur_len;
    write_read(new_rec,in_al_hdr,outSAM,outSAM_header);
//    if(std::strcmp(bam_get_qname(in_rec),"AF004885_119_269_0:0:2_0:0:2_78/1")==0) {
//        std::cout<<bam_get_qname(in_rec)<<std::endl;
//        std::cout<<cur_len<<std::endl;
//        std::cout<<cur_local_start<<std::endl;
//
//        print_seq(in_rec);
//        print_seq(new_rec);
//        print_cigar(in_rec);
//        print_cigar(new_rec);
//    }
}

void MSA::realign(std::string in_sam,std::string out_sam){
    samFile *msa_hdr_fp = hts_open(this->msa_header_fname.c_str(),"r");
    bam_hdr_t *msa_hdr = sam_hdr_read(msa_hdr_fp);

    samFile *outSAM=sam_open(out_sam.c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_init();
    outSAM_header=bam_hdr_dup(msa_hdr);
    sam_hdr_write(outSAM,outSAM_header);
    bam_hdr_destroy(outSAM_header);

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
            write_read(in_rec,in_al_hdr,outSAM,outSAM_header);
        }
    }

    bam_destroy1(in_rec);
    sam_close(in_al);
    sam_close(outSAM);
}

