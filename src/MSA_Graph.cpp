//
// Created by sparrow on 5/2/19.
//

#include "MSA_Graph.h"

MSA_Graph::MSA_Graph() {
    this->IUPAC = std::unordered_map<std::string,std::string>({{"A","A"},
                                                               {"C","C"},
                                                               {"G","G"},
                                                               {"T","T"},
                                                               {"N","N"},
                                                               {"AG","R"},{"GA","R"},
                                                               {"CT","Y"},{"TC","Y"},
                                                               {"CG","S"},{"GC","S"},
                                                               {"AT","W"},{"TA","W"},
                                                               {"GT","K"},{"TG","K"},
                                                               {"AC","M"},{"CA","M"},
                                                               {"CGT","B"},{"CTG","B"},{"GTC","B"},{"GCT","B"},{"TCG","B"},{"TGC","B"},
                                                               {"AGT","D"},{"ATG","D"},{"GAT","D"},{"GTA","D"},{"TAG","D"},{"TGA","D"},
                                                               {"ACT","H"},{"ATC","H"},{"CAT","H"},{"CTA","H"},{"TCA","H"},{"TAC","H"},
                                                               {"ACG","V"},{"AGC","V"},{"CAG","V"},{"CGA","V"},{"GCA","V"},{"GAC","V"},
                                                               {"ACGT","N"},{"ACTG","N"},{"AGCT","N"},{"AGTC","N"},{"ATCG","N"},{"ATGC","N"},
                                                               {"CAGT","N"},{"CATG","N"},{"CGAT","N"},{"CGTA","N"},{"CTAG","N"},{"CTGA","N"},
                                                               {"GACT","N"},{"GATC","N"},{"GCAT","N"},{"GCTA","N"},{"GTAC","N"},{"GTCA","N"},
                                                               {"TACG","N"},{"TAGC","N"},{"TCAG","N"},{"TCGA","N"},{"TGAC","N"},{"TGCA","N"}});

    this->IUPAC_REV = std::unordered_map<std::string,std::string>({{"A","A"},{"a","A"},
                                                                   {"C","C"},{"c","C"},
                                                                   {"G","G"},{"g","G"},
                                                                   {"T","T"},{"t","T"},
                                                                   {"R","AG"},{"r","AG"},
                                                                   {"Y","CT"},{"y","CT"},
                                                                   {"S","CG"},{"s","CG"},
                                                                   {"W","AT"},{"w","AT"},
                                                                   {"K","GT"},{"k","GT"},
                                                                   {"M","AC"},{"m","AC"},
                                                                   {"B","CGT"},{"b","CGT"},
                                                                   {"D","AGT"},{"d","AGT"},
                                                                   {"H","ACT"},{"h","ACT"},
                                                                   {"V","ACG"},{"v","ACG"},
                                                                   {"N","ACGT"},{"n","ACTG"}});
}

MSA_Graph::MSA_Graph(int length,int num_refs) {
    this->IUPAC = std::unordered_map<std::string,std::string>({{"A","A"},
                                                               {"C","C"},
                                                               {"G","G"},
                                                               {"T","T"},
                                                               {"N","N"},
                                                               {"AG","R"},{"GA","R"},
                                                               {"CT","Y"},{"TC","Y"},
                                                               {"CG","S"},{"GC","S"},
                                                               {"AT","W"},{"TA","W"},
                                                               {"GT","K"},{"TG","K"},
                                                               {"AC","M"},{"CA","M"},
                                                               {"CGT","B"},{"CTG","B"},{"GTC","B"},{"GCT","B"},{"TCG","B"},{"TGC","B"},
                                                               {"AGT","D"},{"ATG","D"},{"GAT","D"},{"GTA","D"},{"TAG","D"},{"TGA","D"},
                                                               {"ACT","H"},{"ATC","H"},{"CAT","H"},{"CTA","H"},{"TCA","H"},{"TAC","H"},
                                                               {"ACG","V"},{"AGC","V"},{"CAG","V"},{"CGA","V"},{"GCA","V"},{"GAC","V"},
                                                               {"ACGT","N"},{"ACTG","N"},{"AGCT","N"},{"AGTC","N"},{"ATCG","N"},{"ATGC","N"},
                                                               {"CAGT","N"},{"CATG","N"},{"CGAT","N"},{"CGTA","N"},{"CTAG","N"},{"CTGA","N"},
                                                               {"GACT","N"},{"GATC","N"},{"GCAT","N"},{"GCTA","N"},{"GTAC","N"},{"GTCA","N"},
                                                               {"TACG","N"},{"TAGC","N"},{"TCAG","N"},{"TCGA","N"},{"TGAC","N"},{"TGCA","N"}});

    this->IUPAC_REV = std::unordered_map<std::string,std::string>({{"A","A"},{"a","A"},
                                                                   {"C","C"},{"c","C"},
                                                                   {"G","G"},{"g","G"},
                                                                   {"T","T"},{"t","T"},
                                                                   {"R","AG"},{"r","AG"},
                                                                   {"Y","CT"},{"y","CT"},
                                                                   {"S","CG"},{"s","CG"},
                                                                   {"W","AT"},{"w","AT"},
                                                                   {"K","GT"},{"k","GT"},
                                                                   {"M","AC"},{"m","AC"},
                                                                   {"B","CGT"},{"b","CGT"},
                                                                   {"D","AGT"},{"d","AGT"},
                                                                   {"H","ACT"},{"h","ACT"},
                                                                   {"V","ACG"},{"v","ACG"},
                                                                   {"N","ACGT"},{"n","ACTG"}});

    this->length = length;
    this->num_refs = num_refs;

    // initialize vertices
    for(int i=0;i<length;i++){
        vertices.insert(MSA_Vertex(num_refs,i));
    }

    // initialize the vector of removed
    this->removed = std::vector<int>(length,0);
    this->removed_cleanup = std::vector<int>(length,0);
}

// this function adds a reference name to the index and create a unique ID
uint16_t MSA_Graph::add_ref(std::string ref_name) {
    return this->index.add_ref(ref_name);
}

void MSA_Graph::add_ref(std::string ref_name, int ref_id){
    this->index.add_ref(ref_name,ref_id);
}

// this function sets a mapping between old and new positions in the refence/MSA
void MSA_Graph::add_pos(uint16_t id, uint32_t new_pos) {
    this->index.add(id,new_pos);
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
    int cl = 0;
    for(int i=0;i<this->vertices.size();i++){
        if(this->removed[i]==0 && this->removed_cleanup[i]==0){
            std::string nt_str = "";
            mv = this->vertices.get(i);
            mv->get_supported_nt_string(nt_str);
            iupac_nt = this->IUPAC[nt_str];
            merged_fp<<iupac_nt;
            nt_str.clear();
            if((cl+1) % 60 == 0){
                merged_fp<<std::endl;
            }
            cl++;
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

void MSA_Graph::find_location(int refID, int ref_start, int end, int& new_start, int& s){
    int sum_removed = std::accumulate(this->removed.begin(), this->removed.begin()+ref_start, 0);
    if(!this->removed[ref_start]){
        new_start = ref_start - sum_removed;
        return;
    }
    else if(sum_removed == this->length){
        new_start = 0;
        s = 0;
        return;
    }
    else{
        new_start = ref_start-sum_removed;
        int next_pos = ref_start;
        s=0;
        while(next_pos < end){
            next_pos = this->vertices.get(next_pos)->get_next_pos4ref(refID);
            s++;
            if(!removed[next_pos]){
                sum_removed = std::accumulate(this->removed.begin(),this->removed.begin()+next_pos,0);
                new_start =  next_pos - sum_removed;
                return;
            }
        }
        new_start = ref_start - sum_removed;
    }
}

int MSA_Graph::get_num_clean_removed(int pos){
    return std::accumulate(this->removed_cleanup.begin(),this->removed_cleanup.begin()+pos,0);
}

void MSA_Graph::fit_read(int refID,int ref_start,int end,int& new_start, int& s, std::vector<int>& not_removed, std::unordered_set<int>& added){ // the last four parameters are the return
    s=0;
    this->find_location(refID,ref_start,end,new_start,s);
    if(s>0){
        for(int i=0;i<s;i++){
            added.insert(i);
        }
    }
    std::vector<int> to_remove;
    int pos_tracker = 0,next_vID=ref_start,cur_vID;
    MSA_Vertex* v,next_v;

    while(true){
        cur_vID = next_vID;
        v = this->vertices.get(cur_vID);
        next_vID = v->get_next_pos4ref(refID);
        if(this->removed[cur_vID]==1){ // current position has been removed - causes insertion
            added.insert(pos_tracker);
        }
        else{
            v->inc_ref(refID);
            v->set_mapped();
        }
        int tmp_pos = pos_tracker;
        for(int i=cur_vID+1;i<next_vID;i++){
            if(this->removed[i]==0 && this->vertices.get(i)->is_mapped()){ // can not be removed - causes deletion
                not_removed.emplace_back(tmp_pos);
                tmp_pos++;
            }
            else{
                to_remove.emplace_back(i);
            }
        }

        if(next_vID>=end){
            break;
        }
        pos_tracker++;
    }
    if(to_remove.size()>0){
        for(auto& vid : to_remove){
            this->removed[vid]=1;
        }
    }
}

int MSA_Graph::get_gff_pos(int refID,int pos){
    int ref_start = get_new_position(refID,pos);

    int sum_removed = std::accumulate(this->removed.begin(), this->removed.begin()+ref_start, 0);
    if(removed[ref_start]){
        std::cerr<<"position has been removed but is still being reported"<<std::endl;
    }
    return ref_start-sum_removed;
}

int MSA_Graph::get_first_pos(int refID){
    return this->index.get_first_pos(refID);
}

int MSA_Graph::get_last_pos(int refID){
    return this->index.get_last_pos(refID);
}

void MSA_Graph::set_removed_cleanup(int start, int end){
    for(int i=start;i<end;i++){
        this->removed_cleanup[i]=1;
    }
}

void MSA_Graph::set_removed_cleanup(int pos){
    this->removed_cleanup[pos]=1;
}

int MSA_Graph::get_most_abundant_refID(int pos,int&refID){
    return this->vertices.get(pos)->get_most_abundant_refID(refID);
}

void MSA_Graph::get_first_mapped_pos(int& pos,int& refID){
    for(int i=0;i<this->vertices.size();i++){
        if(this->vertices.get(i)->is_mapped()){
            // get the most abundant reference at the position. if two equally good exist - return the first one
            pos = i;
            get_most_abundant_refID(pos,refID);
            return;
        }
    }
}

void MSA_Graph::get_last_mapped_pos(int& pos,int& refID){
    for(int i=this->vertices.size()-1;i>0;i--){
        if(this->vertices.get(i)->is_mapped()){
            // get the most abundant reference at the position. if two equally good exist - return the first one
            pos = i;
            get_most_abundant_refID(pos,refID);
            return;
        }
    }
}

void MSA_Graph::clean_gaps(int start,int end) { // find gaps and remove any vertices that do not make sense for consensus
    int gap_start = 0;
    int gap_end = 0;
    MSA_Vertex* mv;
    for(int i=start;i<end;i++){
        if(this->removed[i]){
            continue;
        }
        if(gap_start==0){ // gap not found yet
            if(!this->vertices.get(i)->is_mapped()){
                gap_start = i;
            }
        }
        else if(gap_start && !gap_end){ // gap found but need to search for the end
            if(this->vertices.get(i)->is_mapped()){
                gap_end = i-1;
            }
        }
        else{ // both ends are found - can now evaluate and decide what to do
            std::map<int,int> best_refs;
            std::pair<std::map<int,int>::iterator,bool> br_it;
            int refID;
            for(int j=gap_start;j>(gap_start-10)&&j>0;j--){ // evaluate 10 nearest bases for the most common reference
                int count = get_most_abundant_refID(j,refID);
                br_it = best_refs.insert(std::make_pair(refID,0));
                br_it.first->second+=count;
            }

            for(int j=gap_end;j>(gap_end+10)&&j<end;j++){ // evaluate 10 nearest bases for the most common reference
                int count = get_most_abundant_refID(j,refID);
                br_it = best_refs.insert(std::make_pair(refID,0));
                br_it.first->second+=count;
            }

            // now find which reference to preserve
            int best_ref = 0;
            int best_count = 0;
            for(auto& v : best_refs){
                if(v.second>best_count){
                    best_ref = v.first;
                    best_count = v.second;
                }
            }

            // lastly modify the paths to remove unwanted elements
            for(int j=gap_start;j<gap_end;j++){
                if(this->vertices.get(j)->has_ref(best_ref)){
                    continue;
                }
                else{
                    this->removed_cleanup[j]=1;
                    mv=this->get_vertex(j);
                    mv->set_mapped();
                    mv->inc_ref(best_ref);
                }
            }
        }
    }
}

// tag notes on the graph that correspond to the entries in the annotation
int MSA_Graph::pre_fit_annotation(std::string in_gff_fname){
    std::ifstream in_gff_fp;
    in_gff_fp.open(in_gff_fname.c_str(),std::ios::in);

    if(!in_gff_fp.good()){
        std::cerr<<"@ERROR::::Cannot open GFF to read"<<std::endl;
        exit(-1);
    }

    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss("");
    std::string ref_name,track,feature,start_s,end_s,score,strand,phase,attrs;
    int start,end;

    while(std::getline(in_gff_fp,line)) { // iterate over references
        if(line.front() == '#'){
            continue;
        }
        ss.str(line);
        ss.clear();

        std::getline(ss, ref_name, '\t');
        int refID = this->index.getID(ref_name);
        if(refID==-1){
            std::cerr<<"@ERROR::::GFF Reference sequence not found in the index"<<std::endl;
            return -1;
        }
        std::getline(ss, track, '\t');
        std::getline(ss, feature, '\t');
        std::getline(ss, start_s, '\t');
        std::getline(ss, end_s, '\t');
        std::getline(ss, score, '\t');
        std::getline(ss, strand, '\t');
        std::getline(ss, phase, '\t');
        std::getline(ss, attrs, '\t');

        int new_start = this->get_gff_pos(refID,std::atoi(start_s.c_str())-1);
        int new_end = this->get_gff_pos(refID,std::atoi(end_s.c_str())-1);

        MSA_Vertex* ms = this->vertices.get(new_start);
        ms->inc_ref(refID);
        ms->set_mapped();

        MSA_Vertex* me = this->vertices.get(new_end);
        me->inc_ref(refID);
        me->set_mapped();
    }
    in_gff_fp.close();

    std::cerr << "@LOG::loaded the annotation"<<std::endl;
    return 1;
}

int MSA_Graph::fit_annotation(std::string in_gff_fname,std::string out_gff_fname){
    std::ifstream in_gff_fp;
    in_gff_fp.open(in_gff_fname.c_str(),std::ios::in);

    if(!in_gff_fp.good()){
        std::cerr<<"@ERROR::::Cannot open GFF to read"<<std::endl;
        exit(-1);
    }

    std::ios::sync_with_stdio(false);
    std::string line;
    std::stringstream ss("");
    std::string ref_name,track,feature,start_s,end_s,score,strand,phase,attrs;
    int start,end;

    std::ofstream out_gff_fp(out_gff_fname.c_str());

    while(std::getline(in_gff_fp,line)) { // iterate over references
        if(line.front() == '#'){
            out_gff_fp<<line<<std::endl;
            continue;
        }
        ss.str(line);
        ss.clear();

        std::getline(ss, ref_name, '\t');
        int refID = this->index.getID(ref_name);
        if(refID==-1){
            std::cerr<<"@ERROR::::GFF Reference sequence not found in the index"<<std::endl;
            return -1;
        }
        std::getline(ss, track, '\t');
        std::getline(ss, feature, '\t');
        std::getline(ss, start_s, '\t');
        std::getline(ss, end_s, '\t');
        std::getline(ss, score, '\t');
        std::getline(ss, strand, '\t');
        std::getline(ss, phase, '\t');
        std::getline(ss, attrs, '\t');

        int new_start = this->get_gff_pos(refID,std::atoi(start_s.c_str())-1);
        int new_end = this->get_gff_pos(refID,std::atoi(end_s.c_str())-1);

        out_gff_fp<<"MSA"<<"\t"
                  <<track<<"\t"
                  <<feature<<"\t"
                  <<new_start<<"\t"
                  <<new_end<<"\t"
                  <<score<<"\t"
                  <<strand<<"\t"
                  <<phase<<"\t"
                  <<attrs<<std::endl;

    }

    out_gff_fp.close();
    in_gff_fp.close();

    std::cerr << "@LOG::loaded the annotation"<<std::endl;
    return 1;
}

void MSA_Graph::init_refcouts(){
    if(this->length==0){
        std::cerr<<"empty graph"<<std::endl;
        std::exit(-1);
    }
}

// clean the graph for a projection
void MSA_Graph::set_used(int refid){
    // iterate over all the nodes
    MSA_Vertex* mv;
    bool ref_exists;
    for(int i=0;i<this->length;i++){
        mv=this->get_vertex(i);
        if(mv->has_ref(refid)){
            mv->set_mapped();
            mv->inc_ref(refid);
        }
        else{
            this->removed[i]=1;
        }
    }
}

void MSA_Graph::set_used(int start, int end, int refid){ // for each base in this increment, sets the base as used
    MSA_Vertex* mv;
    for(int i=start;i<end;i++){
        mv=this->get_vertex(i);
        mv->set_mapped();
        mv->inc_ref(refid);
    }
}