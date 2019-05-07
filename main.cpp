#include <iostream>
#include <string>
#include <cstring>

#include "src/arg_parse.h"
#include "src/MSA.h"
#include "src/MGraph.h"

void print_help(){
    std::cout<<"help page"<<std::endl;
}

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

    return 0;
}

int vest_build(int argc,char* argv[]){
    enum Opt_Build {MUS_FP   = 'i',
                    MUS_DB   = 'o'};

    ArgParse args_build("vest_build");
    args_build.add_string(Opt_Build::MUS_FP,"mus","","");
    args_build.add_string(Opt_Build::MUS_DB,"out","","");

    args_build.parse_args(argc,argv);

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
    else if (strcmp(argv[1],"help") == 0 || strcmp(argv[1],"--help") == 0){
        print_help();
    }
    else{
        std::cout<<"Unrecognized Mode: please consult the help manual"<<std::endl;
        print_help();
    }

    return 0;
}