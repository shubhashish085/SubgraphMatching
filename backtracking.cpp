//
// Created by antu on 2/8/24.
//

#include "backtracking.h"
#include <map>
#include <queue>

bool
AlgorithmStore :: DFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    ui neighbor_count = 0;
    VertexID * neighbors =  graph->getVertexNeighbors(vertex_to_be_explored, neighbor_count);

    for(ui iterator = 0; iterator < neighbor_count; iterator++){
        if(!visited[neighbors[iterator]]){
            visited[neighbors[iterator]] = true;
            matching_order.push_back(neighbors[iterator]);
            DFSTraversal(graph, neighbors[iterator], non_tree_edges, matching_order, visited);
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

    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    std::queue<VertexID> bfs_queue;
    bfs_queue.push(start_node);

    while(!bfs_queue.empty()){

        VertexID front_vertex = bfs_queue.pop();
        VertexID * neighbors = graph ->getVertexNeighbors(front_vertex, neighbor_count);
        ui neighbor_count = 0;
        for(ui iterator = 0; iterator < neighbor_count; iterator++){
            if(!visited[neighbors[iterator]]){
                bfs_queue.push(neighbors[iterator]);
                matching_order.push_back(neighbors[iterator]);
                visited[neighbors[iterator]] = true;
            }else{
                pair<VertexID, VertexID> nte = make_pair(front_vertex, neighbors[iterator]);
                non_tree_edges.push_back(nte);
            }
        }
    }

    return true;
}


