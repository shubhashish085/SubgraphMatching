//
// Created by antu on 2/8/24.
//

#include "backtracking.h"
#include <map>


bool
AlgorithmStore :: DFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    ui neighbor_count = 0;
    VertexID * neighbors =  graph->getVertexNeighbors(vertex_to_be_explored, neighbor_count);

    for(ui iterator = 0; iterator < neighbor_count; iterator++){
        if(!visited[iterator]){
            matching_order.push_back(neighbors[iterator]);
            DFSTraversal(graph, iterator, non_tree_edges, matching_order, visited);
        } else{
            pair<VertexID, VertexID> nte = make_pair(vertex_to_be_explored, neighbors[iterator]);
            non_tree_edges.push_back(nte);
        }
    }

    return true;
}

bool
AlgorithmStore :: BFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    return true;
}


