

#include "PruningConstraints.h"
#include <tuple>
#define INVALID_VERTEX_ID 100000000

void PruningConstraints::buildCandidateMap(const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order, std::map<VertexID, std::vector<VertexID>>& candidate_map){

    for(ui i = 0; i < query_graph->getVerticesCount(); i++){
        ui vertex = order[i];
        for(ui j = 0; j < candidates_count[vertex]; j++){
            std::map<VertexID, std::vector<VertexID>>::iterator itr = candidate_map.find(candidates[vertex][j]);

            if(itr != candidate_map.end()){
                (itr->second).push_back(vertex);

            }else{
                std::vector<VertexID> mapping_vtx;

                mapping_vtx.push_back(vertex);
                candidate_map[candidates[vertex][j]] = mapping_vtx;
            }
        }
    }
}


bool PruningConstraints::isCandidate(VertexID query_vertex, VertexID vertexToBeChecked, std::map<VertexID, std::vector<VertexID>>& candidate_map){

    if(candidate_map.find(vertexToBeChecked) == candidate_map.end()){
        return false;
    }

    std::map<VertexID, std::vector<VertexID>>::iterator itr = candidate_map.find(vertexToBeChecked);

    for(auto vtr_element: itr->second){
        if(vtr_element == query_vertex){
            return true;
        }
    }

    return false;
}



void PruningConstraints::checkingCycle(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                        std::vector<std::tuple<VertexID, VertexID>>& edges, std::map<VertexID, std::vector<VertexID>>& candidate_map){
    
    std::cout << "Checking Cycles ..." << std::endl;

    std::vector<std::tuple<VertexID, VertexID, short>> previous_level_tuples;
    std::vector<std::tuple<VertexID, VertexID, short>> current_level_tuples;

    ui* neighbors;
    VertexID u, v;
    ui neighbor_count, pruned_vertex_count = 0;
    bool cycle_found = false;

    for(ui i = 0; i < candidates_count[std::get<0>(edges[0])]; i++){

        cycle_found = false;

        for(ui s = 0; s < edges.size(); s++){
            u = std::get<0>(edges[s]);
            v = std::get<1>(edges[s]);

            neighbor_count = 0;

            if(s == 0){

                neighbors = data_graph->getVertexNeighbors(candidates[u][i], neighbor_count);

                for(ui j = 0; j < neighbor_count; j++){
                    if(isCandidate(v, neighbors[j], candidate_map)){
                        current_level_tuples.push_back(std::make_tuple(neighbors[j], candidates[std::get<0>(edges[0])][i], s + 1));
                    }
                }

                previous_level_tuples.clear();
                previous_level_tuples.assign(current_level_tuples.begin(), current_level_tuples.end());
                current_level_tuples.clear();

            }else{

                for(ui k = 0; k < previous_level_tuples.size(); k++){
                    neighbors = data_graph->getVertexNeighbors(std::get<0>(previous_level_tuples[k]), neighbor_count);

                    for(ui j = 0; j < neighbor_count; j++){
                        if(isCandidate(v, neighbors[j], candidate_map)){
                            current_level_tuples.push_back(std::make_tuple(neighbors[j], candidates[std::get<0>(edges[0])][i], s + 1));
                            if(s + 1 == edges.size() && neighbors[j] == candidates[std::get<0>(edges[0])][i]){
                                cycle_found = true;
                                break;
                            }
                        }
                    }

                    if(cycle_found){
                        break;
                    }
                }

                if(cycle_found){
                    break;
                }

                previous_level_tuples.clear();
                previous_level_tuples.assign(current_level_tuples.begin(), current_level_tuples.end());
                current_level_tuples.clear();
            }
            
        }

        if(!cycle_found){
            candidates[std::get<0>(edges[0])][i] = INVALID_VERTEX_ID;
            pruned_vertex_count++;
        }

        previous_level_tuples.clear();
        current_level_tuples.clear();
    }

    ui next_position = 0;

    for (ui i = 0; i < candidates_count[std::get<0>(edges[0])]; ++i) {

        if (candidates[std::get<0>(edges[0])][i] != INVALID_VERTEX_ID){
            candidates[std::get<0>(edges[0])][next_position++] = candidates[std::get<0>(edges[0])][i];
        }           
    }

    candidates_count[std::get<0>(edges[0])] = next_position; 

    std::cout << "Pruned Vertex Count in Cycle Checking Phase : " << pruned_vertex_count << std::endl;

}