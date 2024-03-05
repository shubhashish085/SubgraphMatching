//
// Created by antu on 2/12/24.
//

#include "graph.h"
#include "backtracking.h"

bool filter_by_neighborhood_label_count(std::unordered_map<LabelID, ui>& d_vtx_nlc, std::unordered_map<LabelID, ui>& q_vtx_nlc){

    //std::cout << "------------------------ In Neighborhood Label Filter Method  ------------------------"  << std::endl;

    if(q_vtx_nlc.size() > d_vtx_nlc.size()){
        return false;
    }

    for (auto iterator : q_vtx_nlc) {
        if(q_vtx_nlc[iterator.first] > d_vtx_nlc[iterator.first]){
            return false;
        }
    }

    return true;
}

bool checkNonTreeEdges(Graph* data_graph, ui* query_match_order_idx, std::vector<std::pair<VertexID, VertexID>>& query_ntes,  std::vector<VertexID>& partial_result){

    //std::cout << "------------------------ In Non Tree Edge Method ------------------------"  << std::endl;

    bool edgeExists = false;

    for(ui i = 0; i < query_ntes.size(); i++){
         std::pair<VertexID, VertexID> nte = query_ntes[i];
         edgeExists = data_graph -> checkEdgeExistence(partial_result[nte.first], partial_result[nte.second]);
         if(!edgeExists){
             return false;
         }
    }

    return true;
}

bool edgeTaken(VertexID start, VertexID end, std::vector<std::pair<VertexID, VertexID>>& all_edges){

    for(ui i = 0; i < all_edges.size(); i++){
        if((all_edges[i].first == start && all_edges[i].second == end) || (all_edges[i].first == end && all_edges[i].second == start)){
            return true;
        }
    }

    return false;
}

void printMatch(int query_vertex_count, ui* c_size, ui* iter, VertexID** c){

    std::cout << "Matched : ";

    for(ui i = 0; i < query_vertex_count; i++){
        std::cout << c[i][iter[i]] << " ";
    }

    std::cout << std::endl;
}

void vectorCopyToArray(std::vector<VertexID>& candidate_vtr, VertexID* vertex_list){

    for(ui i = 0; i < candidate_vtr.size(); i++){
        vertex_list[i] = candidate_vtr[i];
    }
}

void stackBasedDFS(Graph* data_graph, Graph* query_graph, std::vector<ui>& matching_order, std::vector<VertexID>* candidate_vtx_vector){

    std::cout << "######### In Stack Based DFS Method #############" << std::endl;

    ui* c_size = new ui[query_graph -> getVerticesCount()];
    ui* iter = new ui[query_graph -> getVerticesCount()];
    VertexID** c = new VertexID*[query_graph -> getVerticesCount()];

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        c[i] = new VertexID[data_graph-> getGraphMaxDegree()];
        c_size[i] = 0;
    }

    ui l = 0;

    //int count = 0;
    while (true){
        if(l < query_graph -> getVerticesCount()){

            if(l == 0 && iter[l] == c_size[l] && iter[l] != 0){
                break;
            }

            if(c_size[l] == 0){
                vectorCopyToArray(candidate_vtx_vector[matching_order[l]], c[l]);
                //std::cout << "------------- Copied : " << l << std::endl;

                if(candidate_vtx_vector[matching_order[l]].size() > 0){
                    c_size[l] = candidate_vtx_vector[matching_order[l]].size();
                } else{
                    c_size[l] = 0;
                }
                if(l == 0 && c_size[l] == 0) break;
                iter[l] = 0;
            }

            //std::cout << " l : " << l << " iter : " << iter[l] << std::endl;
            if(iter[l] < (c_size[l])){
                l++;
            } else{
                c_size[l] = 0;
                if(l > 0){
                    l--;
                    iter[l]++;
                    //std::cout << "Iter incremented" << std::endl;
                }
            }

        }else{
            printMatch(query_graph -> getVerticesCount(), c_size, iter, c);
            l--;
            iter[l]++;
        }

        /*for(ui j = 0; j < l; j++){

            std::cout << "Iter : " << iter[j] << "C : ";
            for(int k = 0; k < 3; k++){
                std::cout << c[j][k] << " ";
            }
            std::cout << " C size :" << c_size[j] <<std::endl;

        }*/



        std::cout << "--------------------------" << std::endl;

    }
}


void match(Graph* data_graph, Graph* query_graph, VertexID candidate_vtx, std::vector<ui>& matching_order, ui matching_idx,
         std::vector<std::pair<VertexID, VertexID>>& query_ntes,  std::vector<VertexID>& partial_result, std::vector<std::pair<VertexID, VertexID>>& all_edges){

    if(partial_result.size() == query_graph->getVerticesCount()){

        /*bool nte_match = checkNonTreeEdges(data_graph, query_graph->getMatchingOrderIndex(), query_ntes, partial_result);

        if(!nte_match){
            return;
        }*/

        std::cout << "Matched" << std::endl;

        for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
            std::cout << partial_result[i] << " ";
        }
        std::cout << std::endl;

        return;

    }else{

        ui neighbor_count;
        VertexID * neighbors = data_graph -> getVertexNeighbors(candidate_vtx, neighbor_count);
        /*std::cout << "Candidate vertex : " << candidate_vtx << " Neighbor Count : " << neighbor_count << std::endl;

        for(ui j = 0; j < neighbor_count; j++) {
            std::cout << " Neighbor : " << neighbors[j] << " ";
        }
        std::cout << std::endl;*/

        for(ui i = 0; i < neighbor_count; i++){

            //bool isEligible = filter_by_neighborhood_label_count(data_graph-> getNeighborhoodLabelCount()[neighbors[i]],
                                    //query_graph-> getNeighborhoodLabelCount()[matching_order[matching_idx + 1]]) && (data_graph->getVertexLabel(neighbors[i]) == query_graph->getVertexLabel(matching_order[matching_idx + 1]));

            bool edge_already_taken = edgeTaken(candidate_vtx, neighbors[i], all_edges);

            if(!edge_already_taken){
                matching_idx++;
                partial_result.push_back(neighbors[i]);
                all_edges.push_back(std::make_pair(candidate_vtx, neighbors[i]));
                match(data_graph, query_graph, neighbors[i], matching_order, matching_idx, query_ntes, partial_result, all_edges);
                all_edges.pop_back();
                partial_result.pop_back();
                matching_idx--;
            }
        }

    }
}




void enumerate(Graph* data_graph, Graph* query_graph, std::vector<std::pair<VertexID, VertexID>>& non_tree_edges, std::vector<ui>& matching_order){

    std::cout << "------------------------ In Enumerate Method ------------------------"  << std::endl;

    VertexID start_vertex = matching_order[0];
    std::vector<std::pair<VertexID, VertexID>> all_edges;
    std::vector<VertexID>* candidate_vtx_vector;
    std::vector<VertexID> partial_result;

    candidate_vtx_vector = new std::vector<VertexID>[query_graph -> getVerticesCount()];

    for(VertexID i = 0; i < query_graph-> getVerticesCount(); i++){
        for(VertexID j = 0; j < data_graph-> getVerticesCount(); j++){
            if(data_graph->getVertexLabel(j) == query_graph->getVertexLabel(i) &&
               filter_by_neighborhood_label_count(data_graph->getNeighborhoodLabelCount()[j], query_graph->getNeighborhoodLabelCount()[i])){
                candidate_vtx_vector[i].push_back(j);
            }
        }
    }

    /*std::cout << "------------------------ Candidate Vertices --------------------------------" << std::endl;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        //std::cout << "Candidate of vertex " << i << " : ";
        for (ui j = 0; j < candidate_vtx_vector[i].size(); ++j) {
            std::cout << " " << candidate_vtx_vector[i][j] ;
        }
        std::cout << std::endl;
    }*/


    /* Stack Based DFS Strategy */

    stackBasedDFS(data_graph, query_graph, matching_order, candidate_vtx_vector);


    /* Recursive Strategy*/

    /*for(ui j = 0; j < candidate_vtx_vector[start_vertex].size(); j++) {
            std::cout << "--------------- Candidate Vertex : " << candidate_vtx_vector[start_vertex][j] << std::endl;
            partial_result.push_back(candidate_vtx_vector[start_vertex][j]);
            match(data_graph, query_graph, candidate_vtx_vector[start_vertex][j], matching_order, 0, non_tree_edges, partial_result, all_edges);
            all_edges.clear();
            partial_result.clear();
    }*/
}

int main(int argc, char** argv) {

    std::string input_query_graph_file = "../basic_query_graph.graph";
    std::string input_data_graph_file = "../basic_data_graph.graph";

    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);

    Graph* data_graph = new Graph();
    data_graph->loadGraphFromFile(input_data_graph_file);

    query_graph->printGraphMetaData();
    data_graph->printGraphMetaData();

    std::vector<ui> matching_order;
    std::vector<std::pair<VertexID, VertexID>> non_tree_edges;

    std::vector<bool> visited;

    std::unordered_map<VertexID, ui>* vertex_map = query_graph -> getNeighborhoodLabelCount();

    std::cout << "Neighborhood Label Count " << std::endl;

    /*for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << "Vertex ID : " << i << std::endl;
        for(auto iterator: query_graph->getNeighborhoodLabelCount()[i]){
            std::cout << "Label ID " << iterator.first << " count : " << query_graph->getNeighborhoodLabelCount()[i][iterator.first] << std::endl;
        }
    }*/


    for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
        visited.push_back(false);
    }

    AlgorithmStore::BFSTraversal(query_graph, 0, non_tree_edges, matching_order, visited);

    query_graph -> setMatchingOrderIndex(matching_order);

    std::cout << "Matching Order : " ;

    for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
        std::cout << matching_order[i] << " ";
    }

    std::cout << std::endl;

    enumerate(data_graph, query_graph, non_tree_edges, matching_order);

}
