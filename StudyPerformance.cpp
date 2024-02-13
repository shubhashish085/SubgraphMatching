//
// Created by antu on 2/12/24.
//

#include "graph.h"
#include "backtracking.h"

bool filter_by_neighborhood_label_count(std::unordered_map<LabelID, ui>& d_vtx_nlc, std::unordered_map<LabelID, ui>& q_vtx_nlc){

    if(q_vtx_nlc.size() >= d_vtx_nlc.size()){
        return false;
    }

    for (auto iterator : q_vtx_nlc) {
        if(q_vtx_nlc.find(iterator.first) >= d_vtx_nlc.find(iterator.first)){
            return false;
        }
    }

    return true;
}

bool checkNonTreeEdges(Graph* data_graph, ui* query_match_order_idx, std::vector<std::pair<VertexID, VertexID>>& query_ntes,  std::vector<VertexID>& partial_result){

    bool edgeExists = false;

    for(ui i = 0; i < query_ntes.size(); i++){
         std::pair<VertexID, VertexID> nte = query_ntes[i];
         edgeExists = data_graph ->checkEdgeExistence(partial_result[query_match_order_idx[nte.first]], partial_result[query_match_order_idx[nte.second]]);
         if(!edgeExists){
             return false;
         }
    }

    return true;
}


void match(Graph* data_graph, Graph* query_graph, VertexID candidate_vtx, std::vector<ui>& matching_order, ui matching_idx,
         std::vector<std::pair<VertexID, VertexID>>& query_ntes,  std::vector<VertexID>& partial_result){

    if(partial_result.size() == query_graph->getVerticesCount()){

        bool nte_match = checkNonTreeEdges(data_graph, query_graph->getMatchingOrderIndex(), query_ntes, partial_result);

        if(!nte_match){
            return;
        }

        std::cout << "Matched" << std::endl;

        for(ui i = 0; i < query_vertices_count; i++){
            std::cout << partial_result[i] << " ";
        }
        std::cout << std::endl;

        return;

    }else{

        ui neighbor_count;
        VertexID * neighbors = data_graph->getVertexNeighbors(candidate_vtx, neighbor_count);

        for(ui i = 0; i < neighbor_count; i++){

            bool isEligible = filter_by_neighborhood_label_count(data_graph-> getNeighborhoodLabelCount()[neighbors[i]],
                                    query_graph-> getNeighborhoodLabelCount()[matching_order[matching_idx + 1]]) && (data_graph->getVertexLabel(neighbors[i]) == query_graph->getVertexLabel(matching_order[matching_idx + 1]));
            if(isEligible){
                matching_idx++;
                partial_result.push_back(neighbors[i]);
                match(data_graph, query_graph, neighbors[i], matching_order, matching_idx, query_ntes, partial_result);
                partial_result.pop_back();
                matching_idx--;
            }
        }

    }
}


void enumerate(Graph* data_graph, Graph* query_graph, std::vector<std::pair<VertexID, VertexID>>& non_tree_edges, std::vector<ui>& matching_order){

    VertexID start_vertex = matching_order[0];
    std::vector<VertexID> candidate_vtx_vector, partial_result;

    for(VertexID i = 0; i < data_graph-> getVerticesCount(); i++){
        if(data_graph->getVertexLabel(i) == query_graph->getVertexLabel(start_vertex) &&
        filter_by_neighborhood_label_count(data_graph->getNeighborhoodLabelCount()[i], query_graph->getNeighborhoodLabelCount()[start_vertex])){
            candidate_vtx_vector.push_back(i);
        }
    }

    for(ui i = 0; i < candidate_vtx_vector.size(); i++){
        match(data_graph, query_graph, candidate_vtx_vector[i], matching_order, 0, non_tree_edges, partial_result);
        partial_result.clear();
    }

}

int main(int argc, char** argv) {

    std::string input_query_graph_file = "query.graph";
    std::string input_data_graph_file = "data.graph";

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);

    Graph* data_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_data_graph_file);

    std::vector<ui> matching_order;
    std::vector<std::pair<VertexID, VertexID>> non_tree_edges;

    std::vector<bool> visited;

    for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
        visited.push_back(false);
    }

    AlgorithmStore::BFSTraversal(query_graph, 0, non_tree_edges, matching_order, visited);

    enumerate(data_graph, query_graph, non_tree_edges, matching_order);

}
