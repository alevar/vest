#include <iostream>
#include <string>
#include <cstring>

#include "src/arg_parse.h"
#include "src/MSA.h"

void print_help(){
    std::cout<<"help page"<<std::endl;
}

// extracts VCF from the Vest-realigned files based on the optional arguments
int vest_vcf(int argc,char* argv[]){
    enum Opt_VCF {IN_SAM = 'i',
        OUT_VCF= 'o'};

    ArgParse args_vcf("vest_inspect");
    args_vcf.add_string(Opt_VCF::IN_SAM,"input","","path to the vest-realigned SAM or BAM file");
    args_vcf.add_string(Opt_VCF::OUT_VCF,"output","","path to the output VCF file");

    args_vcf.parse_args(argc,argv);

    return 0;
}

// allows extraction of the graph-encoded information from the database
int vest_inspect(int argc,char* argv[]){
    enum Opt_Inspect {MSA_DB = 'x',
                      OUT_MSA= 'm'};

    ArgParse args_inspect("vest_inspect");
    args_inspect.add_string(Opt_Inspect::MSA_DB,"db","","path to the vest database");
    args_inspect.add_string(Opt_Inspect::OUT_MSA,"msa","","output file for the MSA encoded in the database");

    args_inspect.parse_args(argc,argv);

    return 0;
}

// performs realignment algorithm
int vest_realign(int argc,char* argv[]){
    enum Opt_Realign {INPUT_FP= 'i',
                    OUTPUT= 'o',
                    MUS_DB = 'x',
                    GFF = 'a'};

    ArgParse args_realign("vest_realign");
    args_realign.add_string(Opt_Realign::INPUT_FP,"input","","path to the input alignment");
    args_realign.add_string(Opt_Realign::OUTPUT,"output","","base name for the output files");
    args_realign.add_string(Opt_Realign::MUS_DB,"db","","database build with vest build from the ultiple sequence alignment");
    args_realign.add_string(Opt_Realign::GFF,"ann","","annotation of genomic features with respect to one of the genomes");

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

    msa.realign(args_realign.get_string(INPUT_FP),args_realign.get_string(OUTPUT));

    return 0;
}

// builds the database
int vest_build(int argc,char* argv[]){
    enum Opt_Build {MUS_FP   = 'i',
                    MUS_DB   = 'o'};

    ArgParse args_build("vest_build");
    args_build.add_string(Opt_Build::MUS_FP,"mus","","");
    args_build.add_string(Opt_Build::MUS_DB,"out","","");

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

    msa.save_graph(args_build.get_string(MUS_DB),cl);

    return 0;
}

int main(int argc, char* argv[]) {

    if(strcmp(argv[1],"build") == 0){
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