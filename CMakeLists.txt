cmake_minimum_required(VERSION 3.14)
project(vest)

set(CMAKE_CXX_STANDARD 11)

add_executable(vest main.cpp src/MSA.cpp src/MSA.h src/MGraph.cpp src/MGraph.h src/arg_parse.cpp src/arg_parse.h src/MSA_List.cpp src/MSA_List.h)