#include <iostream>
#include <string>
#include <cstring>

#include "src/arg_parse.h"
#include "src/MSA.h"

void print_help(){
    std::cerr<<"Vest"<<std::endl;
    std::cerr<<"Modes:"<<std::endl;
    std::cerr<<"\tbuild - build the graph from the MSA"<<std::endl;
    std::cerr<<"\trealign - map alignments onto the MSA and infer consensus"<<std::endl;
    std::cerr<<"\tinspect - inspect the graph"<<std::endl;
}

// extracts VCF from the Vest-realigned files based on the optional arguments
int vest_vcf(int argc,char* argv[]){
    enum Opt_VCF {IN_SAM = 'i',
        OUT_VCF= 'o'};

    ArgParse args_vcf("vest_inspect");
    args_vcf.add_string(Opt_VCF::IN_SAM,"input","","path to the vest-realigned SAM or BAM file",true);
    args_vcf.add_string(Opt_VCF::OUT_VCF,"output","","path to the output VCF file",false);

    args_vcf.parse_args(argc,argv);

    return 0;
}

// allows extraction of the graph-encoded information from the database
int vest_inspect(int argc,char* argv[]){
    enum Opt_Inspect {MSA_DB = 'x',
                      OUT_MSA= 'm'};

    ArgParse args_inspect("vest_inspect");
    args_inspect.add_string(Opt_Inspect::MSA_DB,"db","","path to the vest database",true);
    args_inspect.add_string(Opt_Inspect::OUT_MSA,"msa","","output file for the MSA encoded in the database",false);

    args_inspect.parse_args(argc,argv);

    return 0;
}

// performs realignment algorithm
int vest_realign(int argc,char* argv[]){
    enum Opt_Realign {INPUT_FP= 'i',
                    OUTPUT= 'o',
                    MUS_DB = 'x',
                    GFF = 'a',
                    BED = 'b',
                    GAPFILLNAME = 'n',
                    PROJECTION_NAME = 'p'};

    ArgParse args_realign("vest_realign");
    args_realign.add_string(Opt_Realign::INPUT_FP,"input","","path to the input alignment",false);
    args_realign.add_string(Opt_Realign::OUTPUT,"output","","base name for the output files",true);
    args_realign.add_string(Opt_Realign::MUS_DB,"db","","database build with vest build from the ultiple sequence alignment",true);
    args_realign.add_string(Opt_Realign::GFF,"ann","","annotation of genomic features with respect to one of the genomes",false);
    args_realign.add_string(Opt_Realign::BED,"bed","","bed-formatted interval to be converted to the reference sequence",false);
    args_realign.add_string(Opt_Realign::GAPFILLNAME,"gname","","Name of one of the sequences in the pan genome to be used when resolving a gap in the assembly. If not provided, the default behavior is to use the sequence supported by the reads on both sides of the junction",false);
    args_realign.add_string(Opt_Realign::PROJECTION_NAME,"pname","","Name of the sequence in the pan genome to which the alignments will be projected",false);

    args_realign.parse_args(argc,argv);

    std::string cl="vest ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    MSA msa;
    msa.load_graph(args_realign.get_string(MUS_DB),cl);

    if(args_realign.is_set(Opt_Realign::GAPFILLNAME)){
        msa.set_gapfillname(args_realign.get_string(Opt_Realign::GAPFILLNAME));
    }

    if(args_realign.is_set(Opt_Realign::INPUT_FP)){
        msa.realign(args_realign.get_string(INPUT_FP),args_realign.get_string(OUTPUT));
    }

    if(args_realign.is_set(Opt_Realign::PROJECTION_NAME)){
        msa.set_projection_name(args_realign.get_string(Opt_Realign::GAPFILLNAME));
    }

    if(args_realign.is_set(Opt_Realign::GFF)){
        std::string out_gff_fname = args_realign.get_string(Opt_Realign::OUTPUT);
        out_gff_fname.append(".gff");
        msa.fit_annotation(args_realign.get_string(Opt_Realign::GFF),out_gff_fname);
    }

    if(args_realign.is_set(Opt_Realign::BED)){
        std::string out_bed_fname = args_realign.get_string(Opt_Realign::OUTPUT);
        out_bed_fname.append(".bed");
        msa.fit_bed(args_realign.get_string(Opt_Realign::BED),out_bed_fname);
    }

    return 0;

    // TODO: need to write a script which can extract variants from the assembly/alignment
    // TODO: another use case for this would be to partition viruses based on error profiles into single cells/origins in bulk experiments (similar to phasing)
}

// builds the database
int vest_build(int argc,char* argv[]){
    enum Opt_Build {MUS_FP   = 'i',
                    MUS_DB   = 'o'};

    ArgParse args_build("vest_build");
    args_build.add_string(Opt_Build::MUS_FP,"mus","","path to the multiple sequence alignment",true);
    args_build.add_string(Opt_Build::MUS_DB,"out","","path and basename for the output database",true);

    args_build.parse_args(argc,argv);

    MSA msa(args_build.get_string(MUS_FP));

    std::string cl="vest ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    std::cerr<<"saving the index"<<std::endl;

    msa.save_graph(args_build.get_string(MUS_DB),cl);

    return 0;
}

int main(int argc, char* argv[]) {
    if(argc<=1){
        print_help();
    }
    else if(strcmp(argv[1],"build") == 0){
        std::cerr<<"building index"<<std::endl;
        int argc_build=argc-1;
        char* argv_build[argc_build];
        memcpy(argv_build, argv+1, argc_build*sizeof(char*));
        vest_build(argc_build,argv_build);
    }
    else if(strcmp(argv[1],"realign") ==0 ){
        std::cerr<<"realigning"<<std::endl;
        int argc_realign=argc-1;
        char* argv_realign[argc_realign];
        memcpy(argv_realign, argv+1, argc_realign*sizeof(char*));
        vest_realign(argc_realign,argv_realign);
    }
    else if(strcmp(argv[1],"insepct") ==0 ){
        std::cerr<<"inspecting"<<std::endl;
        int argc_inspect=argc-1;
        char* argv_inspect[argc_inspect];
        memcpy(argv_inspect, argv+1, argc_inspect*sizeof(char*));
        vest_inspect(argc_inspect,argv_inspect);
    }
    else if (strcmp(argv[1],"help") == 0 || strcmp(argv[1],"--help") == 0){
        print_help();
    }
    else{
        std::cout<<"Unrecognized Mode: please consult the help manual"<<std::endl;
        print_help();
    }

    return 0;
}

// TODO: A better way of performing projections
//    1. if we want to project the alignment onto a single genome (eg. K03455) before performing the fitting (during which nodes are being tagged as used), we can pre-tag all nodes on the K03455 path as used, and the rest as unused
//        thus enforcing fitting of the alignments strictly onto the chosen genome
//    2. If we want to use full graph consensus but retain only a single genome in the gaps - we can tag the unused nodes with appropriate label and everything else tag as unused