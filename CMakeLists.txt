cmake_minimum_required(VERSION 3.10)
project(vest)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

set(CMAKE_CXX_STANDARD 11)

include("cmake/zlib.cmake")
include("cmake/htslib.cmake")

add_executable(vest
        main.cpp
        src/MSA.cpp
        src/MSA.h
        src/MSA_Graph.cpp
        src/MSA_Graph.h
        src/arg_parse.cpp
        src/arg_parse.h
        src/MSA_List.cpp
        src/MSA_List.h
        src/MSA_Index.h
        src/MSA_Vertex.h
        src/MSA_Edge.h
        src/FastaTools.cpp
        src/FastaTools.h
        src/codons.h
        src/codons.cpp
        src/GBase.h
        src/GBase.cpp
        src/gdna.h
        src/gdna.cpp
        src/GFaSeqGet.h
        src/GFaSeqGet.cpp
        src/gff.h
        src/gff.cpp
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/tokenize.h
        src/tokenize.cpp)

add_dependencies(vest htslib zlib)
target_link_libraries(vest ${htslib_LIB} ${zlib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl)