//
// Created by antu on 2/12/24.
//

#include "matchingcommand.h"
#include "BuildTable.h"
#include "graph.h"
#include "backtracking.h"
#include "FilterVertices.h"
#include "GeneratingFilterPlan.h"
#include "Enumeration.h"
#include "ParallelEnumeration.h"
#include "wtime.h"
#include <chrono>
#include <limits>
#include <fstream>

#define INVALID_VERTEX_ID 100000000
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))
#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)

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


/*void match(Graph* data_graph, Graph* query_graph, std::vector<ui>& matching_order, std::vector<VertexID>* candidate_vtx_vector){

    if(partial_result.size() == query_graph->getVerticesCount()){

        std::cout << "Matched" << std::endl;

        for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
            std::cout << partial_result[i] << " ";
        }
        std::cout << std::endl;

        return;

    }else{

        ui neighbor_count;
        VertexID * neighbors = data_graph -> getVertexNeighbors(candidate_vtx, neighbor_count);

        for(ui i = 0; i < neighbor_count; i++){

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




void match(Graph* data_graph, Graph* query_graph, VertexID candidate_vtx, std::vector<ui>& matching_order, ui matching_idx,
         std::vector<std::pair<VertexID, VertexID>>& query_ntes,  std::vector<VertexID>& partial_result, std::vector<std::pair<VertexID, VertexID>>& all_edges){

    if(partial_result.size() == query_graph->getVerticesCount()){

        std::cout << "Matched" << std::endl;

        for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
            std::cout << partial_result[i] << " ";
        }
        std::cout << std::endl;

        return;

    }else{

        ui neighbor_count;
        VertexID * neighbors = data_graph -> getVertexNeighbors(candidate_vtx, neighbor_count);

        for(ui i = 0; i < neighbor_count; i++){

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
}*/





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

void exploreByRecursion(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                          ui *candidates_count, ui *order, TreeNode *& tree,
                          size_t output_limit_num, size_t &call_count){

    size_t embedding_count = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];


    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    for(ui i = 0; i < idx_count[cur_depth]; i++){
        Enumerate::exploreRecursiveFashion(data_graph, query_graph, candidates, candidates_count, embedding,
                cur_depth, max_depth, order, idx, idx_count, visited_vertices, valid_candidate,
                tree, embedding_count, call_count);
    }


}


void analyseFiltering(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui filtered_candidate_count = 0;
    bool* candidate_trace = new bool[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);
    std::fill(candidate_trace, candidate_trace + data_graph -> getVerticesCount(), false);

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        for(ui j = 0; j < candidates_count[i]; j++){
            candidate_trace[candidates[i][j]] = true;
        }
    }

    for(ui i = 0; i < data_graph -> getVerticesCount(); i++){
        if(candidate_trace[i]){
            filtered_candidate_count += 1;
        }
    }

    std::cout << " Query Graph : " << std::endl;
    query_graph -> printGraphMetaData();
    std::cout << " Data Graph : " << std::endl;
    data_graph -> printGraphMetaData();
    std::cout << " Filtered Vertex Count : " << filtered_candidate_count << std::endl;

}


void analyseResult(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    int thread_count = 2;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    auto start = std::chrono::high_resolution_clock::now();

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    auto end = std::chrono::high_resolution_clock::now();

    double preprocessing_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    //Stack Based Strategy
    std::cout << "Exploration Started" << std::endl;
    double start_time = wtime();
    embedding_count = Enumerate::exploreAndAnalysis(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, vertex_participating_in_embedding,output_limit, call_count, output_file_path);
    double end_time = wtime();
    std::cout << "Serial Stack Based Strategy Embedding Count : " << embedding_count << " Call Count : " << call_count << std::endl;
    std::cout << "Time " << end_time - start_time << std::endl;

    std::ofstream outputfile;
    outputfile.open(output_file_path, std::ios::app);

    outputfile << "Enumerate time (seconds): " << end_time - start_time << std::endl;
    outputfile << "#Embeddings: " << embedding_count << std::endl;
    outputfile << "Call Count: " << call_count << std::endl;

    outputfile << "-------------------------------- End -----------------------------------" <<std::endl;

    size_t memory = BuildTable::computeMemoryCostInBytes(query_graph, data_graph, candidates_count);
    std::cout << "Amount of Memory Used(MB) : " << BYTESTOMB(memory) << std::endl;
    std::cout << "Preprocessing Time : " << NANOSECTOSEC(preprocessing_time_in_ns) << std::endl;


    outputfile.flush();
    outputfile.close();


}




void analyseDegree(Graph* query_graph, Graph* data_graph){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui loop_count = 4;
    int thread_count[] = {2, 4, 8, 16};


    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    ui* cand_degree_offset = new ui[candidates_count[start_vertex] + 1];
    ui* thread_wise_degree = NULL;
    ui* candidate_limit = NULL;
    ui* allocated_degree_dist = NULL;
    cand_degree_offset[0] = 0;

    ui total_degree = 0;
    ui avg_degree = 0;

    std::cout << "Exploration Started" << std::endl;
    for (ui i = 0; i < loop_count; i++){

        ui avg_size = candidates_count[start_vertex] / thread_count[i];
        //thread_wise_degree = new ui[thread_count[i]];
        candidate_limit = new ui[thread_count[i]];
        //allocated_degree_dist = new ui[thread_count[i]];

        for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
            cand_degree_offset[j] = cand_degree_offset[j - 1] + data_graph->getVertexDegree(candidates[start_vertex][j - 1]);
        }

        total_degree = cand_degree_offset[candidates_count[start_vertex]];
        avg_degree = total_degree / thread_count[i];

        /*for(ui j = 0; j < thread_count[i]; j++){
            if(j == thread_count[i] - 1){
                thread_wise_degree[j] = cand_degree_offset[candidates_count[start_vertex]] - cand_degree_offset[j * avg_size];
            }else {
                thread_wise_degree[j] = cand_degree_offset[(j + 1) * avg_size] - cand_degree_offset[j * avg_size];
            }
        }*/

        ui last_thread_degree_offset = 0;
        ui thread_idx = 0;

        for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
            if(thread_idx == thread_count[i] - 1){
                //allocated_degree_dist[thread_idx] = cand_degree_offset[candidates_count[start_vertex]] - last_thread_degree_offset;
                candidate_limit[thread_idx++] = candidates_count[start_vertex];
                last_thread_degree_offset = cand_degree_offset[candidates_count[start_vertex]];
                break;
            }else if(j == candidates_count[start_vertex]){
                //allocated_degree_dist[thread_idx] = cand_degree_offset[j] - last_thread_degree_offset;
                candidate_limit[thread_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }else if(cand_degree_offset[j] - last_thread_degree_offset >= avg_degree){
                //allocated_degree_dist[thread_idx] = cand_degree_offset[j] - last_thread_degree_offset;
                candidate_limit[thread_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }
        }

        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Number of Threads : " << thread_count[i] << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        for(ui j = 0; j < thread_count[i]; j++){
            std::cout << "Thread No " << j  << " :: even degree summation : " << thread_wise_degree[j] << std::endl;
            std::cout << "Thread No " << j  << " :: candidate limit : " << candidate_limit[j] << " allocated_degree_dist : " << allocated_degree_dist[j] << std::endl;
        }

        std::cout << "-----------------------------------------" << std::endl;

        //delete[] thread_wise_degree;
        //delete[] candidate_limit;
        //delete[] allocated_degree_dist;

    }

}

void analyseParallelizationForWeakScaling(Graph* query_graph, Graph* data_graph, const std::string& output_file_path, int& thread_count){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    VertexID start_vertex = matching_order[0];

    ui* cand_degree_offset = new ui[candidates_count[start_vertex] + 1];
    cand_degree_offset[0] = 0;
    ui* candidate_limit = NULL;

    for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
        cand_degree_offset[j] = cand_degree_offset[j - 1] + data_graph->getVertexDegree(candidates[start_vertex][j - 1]);
    }

    //Parallel Strategy
    double start_time, end_time;

    std::cout << "Exploration Started" << std::endl;


    candidate_limit = new ui[thread_count];

    ui total_degree = cand_degree_offset[candidates_count[start_vertex]];
    ui avg_degree = total_degree / thread_count;

    ui last_thread_degree_offset = 0;
    ui thread_idx = 0;

    for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
        if(thread_idx == thread_count - 1){
            candidate_limit[thread_idx++] = candidates_count[start_vertex];
            last_thread_degree_offset = cand_degree_offset[candidates_count[start_vertex]];
            break;
        }else if(j == candidates_count[start_vertex]){
            candidate_limit[thread_idx++] = j;
            last_thread_degree_offset = cand_degree_offset[j];
        }else if(cand_degree_offset[j] - last_thread_degree_offset >= avg_degree){
            candidate_limit[thread_idx++] = j;
            last_thread_degree_offset = cand_degree_offset[j];
        }
    }

    std::cout << "Candidate Limit : " << std::endl;

    for (ui j = 0; j < thread_idx; j++){
        std::cout << candidate_limit[j] << " " ;
    }
    std::cout << "--------------------------" << std::endl;

    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    size_t** embedding_cnt_array = ParallelEnumeration::exploreWithEvenDegreeDist(data_graph, query_graph, candidates,
                                                                              candidates_count, matching_order, query_tree, output_limit, call_count, thread_count, candidate_limit);
    for(ui idx = 0; idx < thread_count; idx++){
        embedding_count += embedding_cnt_array[idx][0];
    }

    end_time = wtime();
    ParallelEnumeration::writeResult(output_file_path, thread_count, call_count, embedding_count);

    std::cout << "Time " << end_time - start_time << std::endl;


}


void analyseParallelizationWithLoadBalance(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    ui loop_count = 4;
    int thread_count[] = {2, 4, 8, 16};
    //int thread_count[] = {16};
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    VertexID start_vertex = matching_order[0];

    ui* cand_degree_offset = new ui[candidates_count[start_vertex] + 1];
    cand_degree_offset[0] = 0;
    ui* candidate_limit = NULL;

    for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
        cand_degree_offset[j] = cand_degree_offset[j - 1] + data_graph->getVertexDegree(candidates[start_vertex][j - 1]);
    }

    //Parallel Strategy
    double start_time, end_time;

    std::cout << "Exploration Started" << std::endl;
    for (ui i = 0; i < loop_count; i++){

        candidate_limit = new ui[thread_count[i]];

        ui total_degree = cand_degree_offset[candidates_count[start_vertex]];
        ui avg_degree = total_degree / thread_count[i];

        ui last_thread_degree_offset = 0;
        ui thread_idx = 0;

        for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
            if(thread_idx == thread_count[i] - 1){
                candidate_limit[thread_idx++] = candidates_count[start_vertex];
                last_thread_degree_offset = cand_degree_offset[candidates_count[start_vertex]];
                break;
            }else if(j == candidates_count[start_vertex]){
                candidate_limit[thread_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }else if(cand_degree_offset[j] - last_thread_degree_offset >= avg_degree){
                candidate_limit[thread_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }
        }

        std::cout << "Candidate Limit : " << std::endl;

        for (ui j = 0; j < thread_idx; j++){
            std::cout << candidate_limit[j] << " " ;
        }
        std::cout << "--------------------------" << std::endl;

        embedding_count = 0;
        call_count = 0;

        start_time = wtime();
        size_t** embedding_cnt_array = ParallelEnumeration::exploreWithEvenDegreeDist(data_graph, query_graph, candidates,
                                   candidates_count, matching_order, query_tree, output_limit, call_count, thread_count[i], candidate_limit);
        /*size_t** embedding_cnt_array = ParallelEnumeration::exploreWithDynamicLoadBalance(data_graph, query_graph, candidates,
                                                                                      candidates_count, matching_order, query_tree, output_limit, call_count, thread_count[i]);*/
        for(ui idx = 0; idx < thread_count[i]; idx++){
            embedding_count += embedding_cnt_array[idx][0];
        }

        end_time = wtime();
        ParallelEnumeration::writeResult(output_file_path, thread_count[i], call_count, embedding_count);

        std::cout << "Time " << end_time - start_time << std::endl;
    }

}

void analyseParallelizationWithDynamicLoadBalance(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    ui loop_count = 4;
    int thread_count[] = {2, 4, 8, 16};
    //int thread_count[] = {16};
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    //Parallel Strategy
    double start_time, end_time;

    std::cout << "Exploration Started" << std::endl;
    for (ui i = 0; i < loop_count; i++){

        embedding_count = 0;
        call_count = 0;

        start_time = wtime();
        size_t** embedding_cnt_array = ParallelEnumeration::exploreWithDynamicLoadBalance(data_graph, query_graph, candidates,
                                                                                          candidates_count, matching_order, query_tree, output_limit, call_count, thread_count[i]);
        for(ui idx = 0; idx < thread_count[i]; idx++){
            embedding_count += embedding_cnt_array[idx][0];
        }

        end_time = wtime();
        ParallelEnumeration::writeResult(output_file_path, thread_count[i], call_count, embedding_count);

        std::cout << "Time " << end_time - start_time << std::endl;
    }

}

void writeThreadIdAndCount(const std::string& file_path, ui** candidates, ui* candidate_count, ui start_vertex, ui* thread_map, size_t* result){
    std::ofstream outputfile;
    outputfile.open(file_path, std::ios::app);

    for (ui i = 0; i < candidate_count[start_vertex]; i++){
        outputfile << candidates[start_vertex][i] << " "<< thread_map[i] << " " << result[i] << std::endl;

        if(i % 1000 == 0){
            outputfile.flush();
        }

    }

    outputfile.flush();
    outputfile.close();
}

void analyseParallelizationWithDynamicLoadBalanceForAnalysis(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    ui loop_count = 4;
    int thread_count[] = {2, 4, 8, 16};
    //int thread_count[] = {16};
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    //Parallel Strategy
    double start_time, end_time;

    VertexID start_vertex = matching_order[0];

    ui* thread_map = new ui [candidates_count[start_vertex]];
    size_t* result = new size_t [candidates_count[start_vertex]];


    std::cout << "Exploration Started" << std::endl;
    for (ui i = 3; i < loop_count; i++){

        embedding_count = 0;
        call_count = 0;

        start_time = wtime();
        size_t** embedding_cnt_array = ParallelEnumeration::exploreWithDynamicLoadBalanceForAnalysis(data_graph, query_graph, candidates,
                                                                                          candidates_count, matching_order, query_tree, output_limit, call_count, thread_count[i],
                                                                                          thread_map, result);
        for(ui idx = 0; idx < thread_count[i]; idx++){
            embedding_count += embedding_cnt_array[idx][0];
        }

        end_time = wtime();
        writeThreadIdAndCount(output_file_path, candidates, candidates_count, start_vertex, thread_map, result);
        ParallelEnumeration::writeResult(output_file_path, thread_count[i], call_count, embedding_count);

        std::cout << "Time " << end_time - start_time << std::endl;
    }

}


void analyseParallelization(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    ui loop_count = 1;
    //int thread_count[] = {2, 4, 8, 16};
    int thread_count[] = {2};
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    //Parallel Strategy
    double start_time, end_time;

    std::cout << "Exploration Started" << std::endl;
    for (ui i = 0; i < loop_count; i++){

        embedding_count = 0;
        call_count = 0;

        start_time = wtime();
        ui** embedding_cnt_array = ParallelEnumeration::exploreWithPadding(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, output_limit, call_count, thread_count[i]);
        for(ui idx = 0; idx < thread_count[i]; idx++){
            embedding_count += embedding_cnt_array[idx][0];
        }

        end_time = wtime();
        ParallelEnumeration::writeResult(output_file_path, thread_count[i], call_count, embedding_count);

        std::cout << "Time " << end_time - start_time << std::endl;
    }

}




void studyPerformance(Graph* query_graph, Graph* data_graph){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    int thread_count = 2;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;



    //Stack Based Strategy
    std::cout << "Exploration Started" << std::endl;
    double start_time = wtime();
    Enumerate::explore(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, output_limit, call_count);
    double end_time = wtime();
    std::cout << "Serial Stack Based Strategy Embedding Count : " << embedding_count << " Call Count : " << call_count << std::endl;
    std::cout << "Time " << end_time - start_time << std::endl;

    //Recursive Strategy
    //exploreByRecursion(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, output_limit, call_count);
    //std::cout << "Recursive Strategy Embedding Count : " << embedding_count << " Call Count : " << call_count << std::endl;

    //Parallel OpenMP
    /*double start_time = wtime();
    ui* embedding_cnt_array = ParallelEnumeration::explore(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, output_limit, call_count, thread_count);
    double end_time = wtime();
    for(ui i = 0; i < thread_count; i++){
        embedding_count += embedding_cnt_array[i];
    }

    std::cout << "Parallel Strategy Embedding Count : " << embedding_count << std::endl;
    std::cout << "Time" << end_time - start_time << std::endl;*/



    //Recursive Strategy Without Candidate
    /*ui * embedding = new ui[query_graph -> getVerticesCount()];
    bool* visited_vertices = new bool[data_graph -> getVerticesCount()];

    for(ui i = 0; i < data_graph->getVerticesCount(); i++){
        visited_vertices[i] = false;
    }

    start_time = wtime();
    Enumerate::exploreWithoutCandidate(data_graph, query_graph, matching_order, embedding, 0, visited_vertices, query_tree, embedding_count, call_count);
    end_time = wtime();
    std::cout << "Recursive Strategy Without Candidate Embedding Count : " << embedding_count << " Call Count : " << call_count << std::endl;
    std::cout << "Time" << end_time - start_time << std::endl;*/


}


//Main Run
/*
int main(int argc, char** argv) {

    std::string input_query_graph_file = "../tests/basic_query_graph.graph";
    std::string input_data_graph_file = "../tests/basic_data_graph.graph";
    //std::string input_data_graph_file = "../tests/data_graph_4_wo_label.graph";
    //std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/com-dblp.ungraph.txt";
    //std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/soc-LiveJournal1.txt";
    //std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/roadNet-CA.txt";

    Graph* query_graph = new Graph();
    //query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->loadDirectedGraphFromFile(input_query_graph_file);

    query_graph->printGraphMetaData();

    Graph* data_graph = new Graph();
    //data_graph->loadGraphFromFile(input_data_graph_file);
    //data_graph->loadGraphFromFileWithEdge(input_data_graph_file);
    //data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    data_graph->loadDirectedGraphFromFile(input_data_graph_file);

    data_graph->printGraphMetaData();

    std::vector<ui> matching_order;
    std::vector<std::pair<VertexID, VertexID>> non_tree_edges;

    //std::vector<bool> visited;
    //int* parent_vtr;

    std::unordered_map<VertexID, ui>* vertex_map = query_graph -> getNeighborhoodLabelCount();

    std::cout << "Neighborhood Label Count " << std::endl;


    /*for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
        visited.push_back(false);
    }

    studyPerformance(query_graph, data_graph);

}
*/

// Filtering Statistics
int main(int argc, char** argv) {

    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string output_performance_file = command.getOutputFilePath();


    std::cout << " Query Graph : " << input_query_graph_file << std::endl;
    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);
    //query_graph->loadGraphFromFileWithoutStringConversion(input_query_graph_file);

    std::cout << " Data Graph : " << input_data_graph_file << std::endl;
    Graph* data_graph = new Graph();
    data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);

    analyseFiltering(query_graph, data_graph, output_performance_file);
}

/*
int main(int argc, char** argv) {

    /*MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string output_performance_file = command.getOutputFilePath();

    //std::string input_query_graph_file = "../tests/basic_query_graph_wo_label.graph";
    //std::string input_query_graph_file = "../tests/4_node_graph_wo_label.graph";
    //std::string input_query_graph_file = "../tests/5_node_graph_wo_label.graph";
    //std::string input_data_graph_file = "../tests/basic_data_graph_wo_label.graph";
    std::string input_data_graph_file = "/home/antu/Research_Projects/dataset/roadNet-CA.txt";

    //std::string input_query_graph_file = "../tests/basic_query_graph_wo_label.graph";
    //std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/com-lj.ungraph.txt";

    int division_factor[] = {2, 4, 8, 16};

    //std::string output_file = "../analysis/sample_test_parallel_weak_scaling.txt";

    //Graph* query_graph = new Graph();
    //query_graph->loadGraphFromFile(input_query_graph_file);
    //query_graph->loadGraphFromFileWithoutStringConversion(input_query_graph_file);

    Graph* data_graph = new Graph();
    //data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    data_graph->loadDirectedGraphFromFile(input_data_graph_file);
    //analyseParallelizationWithLoadBalance(query_graph, data_graph, output_performance_file);
    //analyseParallelizationWithDynamicLoadBalance(query_graph, data_graph, output_performance_file);
    //analyseResult(query_graph, data_graph, output_performance_file);
    //analyseParallelizationWithDynamicLoadBalanceForAnalysis(query_graph, data_graph, output_performance_file);


    /*for(ui i = 0; i < 4; i++) {

        query_graph->printGraphMetaData();
        data_graph->printGraphMetaData();

        std::vector<ui> matching_order;
        std::vector<std::pair<VertexID, VertexID>> non_tree_edges;

        std::vector<bool> visited;

        for (ui i = 0; i < query_graph->getVerticesCount(); i++) {
            visited.push_back(false);
        }

        //analyseParallelizationForWeakScaling(query_graph, data_graph, output_file, division_factor[i]);

    }

}*/



/*int main(int argc, char** argv) {

    std::string input_query_graph_file = "../tests/basic_query_graph_wo_label.graph";
    //std::string input_query_graph_file = "../tests/4_node_graph_wo_label.graph";
    //std::string input_query_graph_file = "../tests/5_node_graph_wo_label.graph";
    //std::string input_data_graph_file = "../tests/basic_data_graph_wo_label.graph";
    std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/com-lj.ungraph.txt";

    //std::string input_query_graph_file = "../tests/basic_query_graph_wo_label.graph";
    //std::string input_data_graph_file = "../tests/formatted_graph_2048.graph";
    //std::string input_data_graph_file = "../tests/data_graph_4_wo_label.graph";




    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);
    //query_graph->loadGraphFromFileWithoutStringConversion(input_query_graph_file);

    Graph* data_graph = new Graph();
    //data_graph->loadGraphFromFile(input_data_graph_file);
    data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    //data_graph->loadGraphFromFileWithWeight(input_data_graph_file);
    //data_graph->loadDirectedGraphFromFile(input_data_graph_file);

    query_graph->printGraphMetaData();
    data_graph->printGraphMetaData();

    std::vector<ui> matching_order;
    std::vector<std::pair<VertexID, VertexID>> non_tree_edges;

    std::vector<bool> visited;
    int* parent_vtr;

    std::unordered_map<VertexID, ui>* vertex_map = query_graph -> getNeighborhoodLabelCount();

    std::cout << "Neighborhood Label Count " << std::endl;


    for(ui i = 0; i < query_graph-> getVerticesCount(); i++){
        visited.push_back(false);
    }


    analyseResult(query_graph, data_graph, output_file);
    //analyseParallelization(query_graph, data_graph, output_file);
    //analyseParallelizationWithLoadBalance(query_graph, data_graph, output_file);
    //analyseDegree(query_graph, data_graph);


}*/
