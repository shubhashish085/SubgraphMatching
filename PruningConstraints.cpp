

#include "FilterVertices.h"
#include "PruningConstraints.h"
#include <algorithm>
#include <tuple>
#define INVALID_VERTEX_ID 100000000

void PruningConstraints::buildCandidateMap(const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order, std::map<VertexID, std::vector<VertexID>>& candidate_map){

    std::cout << "Building Candidate Map............ " << std::endl;

    candidate_map.clear();
    
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


void PruningConstraints::checkLocalConstraints(const Graph *data_graph, const Graph *query_graph, ui **& candidates, ui *& candidates_count, ui *order, ui step_no){

    std::cout << "Checking Local Constraints ..." << std::endl;

    std::map<VertexID, std::vector<VertexID>>::iterator itr;
    std::map<std::pair<VertexID, VertexID>, bool> vertex_candidacy_map;
    std::map<std::pair<VertexID, VertexID>, bool>::iterator vertex_candidacy_map_itr;

    ui pruned_vtx_count = 0;

    ui* query_neighbors, *data_neighbors;
    ui query_nbr_count = 0, data_nbr_count = 0;
    bool neighbor_match_found = false;

    if(step_no == 0){

        for(ui vtx = 0; vtx < data_graph->getVerticesCount(); vtx++){

            for(ui query_vtx = 0; query_vtx < query_graph -> getVerticesCount(); query_vtx++){

                if(query_graph->getVertexLabel(query_vtx) == data_graph->getVertexLabel(vtx) && query_graph->getVertexDegree(query_vtx) <= data_graph->getVertexDegree(vtx)){
                    vertex_candidacy_map[std::make_pair(vtx, query_vtx)] = true;
                }
            }
        }

    }else{

        for(ui i = 0; i < query_graph->getVerticesCount(); i++){
            for(ui j = 0; j < candidates_count[i]; j++){
                vertex_candidacy_map[std::make_pair(candidates[i][j], i)] = true;
            }
        }

    }
    
    std::cout << "First Stage Selection Done" << std::endl;

    while(true){

        for(ui vtx = 0; vtx < data_graph->getVerticesCount(); vtx++){

            for(ui query_vtx = 0; query_vtx < query_graph->getVerticesCount(); query_vtx++){
                vertex_candidacy_map_itr = vertex_candidacy_map.find(std::make_pair(vtx, query_vtx));

                if(vertex_candidacy_map_itr != vertex_candidacy_map.end()){
                    query_neighbors = query_graph->getVertexNeighbors(query_vtx, query_nbr_count);
                    data_neighbors = data_graph->getVertexNeighbors(vtx, data_nbr_count);

                    neighbor_match_found = false;

                    for(ui i = 0; i < query_nbr_count; i++){
                        for(ui j = 0; j < data_nbr_count; j++){
                            vertex_candidacy_map_itr = vertex_candidacy_map.find(std::make_pair(data_neighbors[j], query_neighbors[i]));

                            if(vertex_candidacy_map_itr != vertex_candidacy_map.end()){
                                neighbor_match_found = true;
                                break;
                            }
                        }

                        if(!neighbor_match_found){
                            vertex_candidacy_map.erase(std::make_pair(vtx, query_vtx));
                            pruned_vtx_count++;
                            break;
                        }
                    }

                    if(!neighbor_match_found){
                        break;
                    }                  
                }
            }
        }
        

        if(pruned_vtx_count == 0){
            break;
        }

        pruned_vtx_count = 0;
    }

    //std::cout << "Pruning Finished" << std::endl;

    ui* tmp_candidate_count = new ui[query_graph->getVerticesCount()];
    std::fill(tmp_candidate_count, tmp_candidate_count + query_graph->getVerticesCount(), 0);


    for(auto it = vertex_candidacy_map.begin(); it != vertex_candidacy_map.end(); ++it){

        candidates[(it->first).second][tmp_candidate_count[(it->first).second]++] = (it->first).first;
    }

    //std::cout << "Candidate Array Build Finished" << std::endl;

    for(ui i = 0; i < query_graph->getVerticesCount(); i++){
        candidates_count[i] = tmp_candidate_count[i];
    }

    //std::cout << "Candidate Count Copy Finished" << std::endl;
}



ui PruningConstraints::checkingCycle(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
                                        std::vector<std::tuple<VertexID, VertexID>>& edges, std::map<VertexID, std::vector<VertexID>>& candidate_map){
    
    std::cout << "Checking Cycles ..." << std::endl;

    std::vector<std::tuple<VertexID, VertexID, short>> previous_level_tuples;
    std::vector<std::tuple<VertexID, VertexID, short>> current_level_tuples;

    std::map<VertexID, std::vector<VertexID>>::iterator candidate_map_itr;

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

            /*candidate_map_itr = candidate_map.find(candidates[std::get<0>(edges[0])][i]);

            if(candidate_map_itr != candidate_map.end()){
                for(ui idx = 0; idx < (candidate_map_itr->second).size(); idx++){
                    if((candidate_map_itr->second)[idx] == std::get<0>(edges[0])){
                        (candidate_map_itr->second).erase((candidate_map_itr->second).begin() + idx);
                    }
                }
            }*/
            
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

    std::cout << "Next Position : " << candidates_count[std::get<0>(edges[0])] << std::endl;
    std::cout << "Pruned Vertex Count in Cycle Checking Phase : " << pruned_vertex_count << std::endl;

    return pruned_vertex_count;

}