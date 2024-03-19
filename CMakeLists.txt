cmake_minimum_required(VERSION 2.8)
project(SubgraphMatchingMain)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3 -g -Wall -march=native -pthread")

set(SUBGRAPH_MATCHING_FILE_LIST
        backtracking.cpp
        backtracking.h
        FilterVertices.h
        FilterVertices.cpp
        GeneratingFilterPlan.h
        GeneratingFilterPlan.cpp
        Enumeration.h
        Enumeration.cpp
        ParallelEnumeration.h
        ParallelEnumeration.cpp
        graph.cpp
        graph.h
        utilities.cpp
        utilities.h
        types.h
        util.h
        wtime.h
        StudyPerformance.cpp)


add_executable(SubgraphMatching.out ${SUBGRAPH_MATCHING_FILE_LIST})

