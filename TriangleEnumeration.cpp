//
// Created by antu on 6/14/24.
//

#include "backtracking.h"
#include "TriangleEnumeration.h"


size_t TriangleEnumeration::enumerateTriangles(const int *data_graph, int output_limit_num, int &call_count) {

    //BFS at all node
    ui vertex_num = data_graph->getVerticesCount();
    int* level = new int[vertex_num];

    for(ui i = 0; i < vertex_num; i++){
        level[i] = -1;
    }

    for(ui i = 0; i < vertex_num; i++){
        if(level[i] == -1){
            AlgorithmStore::bfsTraversal(data_graph, i, level);
        }
    }


    //level wise edge separation and graph building

    //Triangle Counting for same level vertices

    //Triangle counting for different level vertices

}



