//
// Created by antu on 2/8/24.
//

#include "backtracking.h"
#include <map>
#include <queue>

bool
AlgorithmStore :: DFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    std::cout << "------------------------ In DFSTraversal Method ------------------------"  << std::endl;

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
            std::pair<VertexID, VertexID> nte = std::make_pair(vertex_to_be_explored, neighbors[iterator]);
            non_tree_edges.push_back(nte);
        }
    }

    return true;
}

bool
AlgorithmStore :: BFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    std::cout << "------------------------ In BFSTraversal Method ------------------------"  << std::endl;
    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    std::queue<VertexID> bfs_queue;
    bfs_queue.push(start_node);

    std::cout << "Start Node : " << start_node << std::endl;

    std::vector<std::pair<VertexID , VertexID>> visited_edges;

    while(!bfs_queue.empty()){
        ui neighbor_count = 0;
        VertexID front_vertex = bfs_queue.front();
        bfs_queue.pop();
        VertexID * neighbors = graph ->getVertexNeighbors(front_vertex, neighbor_count);

        std::cout << "Neighbor count of " << front_vertex << " : " << neighbor_count << std::endl;

        matching_order.push_back(start_node);

        for(ui iterator = 0; iterator < neighbor_count; iterator++){

            if(!visited[neighbors[iterator]]){
                bfs_queue.push(neighbors[iterator]);
                matching_order.push_back(neighbors[iterator]);
                visited[neighbors[iterator]] = true;
                visited_edges.push_back(std::make_pair(front_vertex, neighbors[iterator]));
            }else{
                bool alreadyExists = false;
                for(ui i = 0; i < visited_edges.size(); i++){
                    if((visited_edges[i].first == front_vertex && visited_edges[i].second == neighbors[iterator]) ||
                    (visited_edges[i].first == neighbors[iterator] && visited_edges[i].second == front_vertex)){
                        alreadyExists = true;
                        break;
                    }
                }

                for(ui i = 0; i < non_tree_edges.size(); i++){
                    if((non_tree_edges[i].first == front_vertex && non_tree_edges[i].second == neighbors[iterator]) ||
                       (non_tree_edges[i].first == neighbors[iterator] && non_tree_edges[i].second == front_vertex)){
                        alreadyExists = true;
                        break;
                    }
                }


                if(!alreadyExists) {

                    std::pair <VertexID, VertexID> nte = std::make_pair(front_vertex, neighbors[iterator]);
                    non_tree_edges.push_back(nte);
                }
            }
        }
    }

    return true;
}


