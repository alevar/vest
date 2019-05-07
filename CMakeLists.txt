cmake_minimum_required(VERSION 3.14)
project(vest)

set(CMAKE_CXX_STANDARD 11)

add_executable(vest main.cpp src/MSA.cpp src/MSA.h src/MSA_Graph.cpp src/MSA_Graph.h src/arg_parse.cpp src/arg_parse.h src/MSA_List.cpp src/MSA_List.h src/MSA_Index.h src/MSA_Vertex.h src/MSA_Edge.h src/FastaTools.cpp src/FastaTools.h)