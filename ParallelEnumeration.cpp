//
// Created by antu on 3/18/24.
//

#include "ParallelEnumeration.h"
#include "Enumeration.h"
#include "wtime.h"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#define NUM_THREADS 8
#define PAD 8


ui** ParallelEnumeration::exploreWithPadding(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                 TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count, int &thread_count){

    std::cout << " ################## explore parallel ##################" << std::endl;

    ui** embedding_cnt_array = new ui* [thread_count];
    double* thread_wise_time = new double [thread_count];

    for(ui i = 0; i < thread_count; i++){
        embedding_cnt_array[i] = new ui[PAD];
    }

    ui* candidate_track = new ui[data_graph->getVerticesCount()];
    ui* candidate_offset = new ui[data_graph->getVerticesCount() + 1];
    ui candidate_csr_count = 0;

    std::fill(candidate_track, candidate_track + data_graph -> getVerticesCount(), 0);

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        candidate_csr_count += candidates_count[i];
    }

    ui* candidate_csr = new ui[candidate_csr_count];

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        for(ui j = 0; j < candidates_count[i]; j++){
            VertexID data_vertex = candidates[i][j];
            candidate_track[data_vertex]++;
        }
    }

    candidate_offset[0] = 0;

    for(ui i = 1; i < data_graph -> getVerticesCount() + 1; i++){
        candidate_offset[i] = candidate_offset[i - 1] + candidate_track[i - 1];
    }


    std::fill(candidate_track, candidate_track + data_graph -> getVerticesCount(), 0);

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        for(ui j = 0; j < candidates_count[i]; j++){
            VertexID data_vertex = candidates[i][j];
            candidate_csr[candidate_offset[data_vertex] + candidate_track[data_vertex]] = i;
            candidate_track[data_vertex]++;
        }
    }

    for (ui i = 0; i < data_graph->getVerticesCount() ; ++i) {
        std::sort(candidate_csr + candidate_offset[i], candidate_csr + candidate_offset[i + 1]); // sorting the query graph parent of every vertex
    }

    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    int avg_size = candidates_count[start_vertex] / thread_count;

    omp_set_num_threads(thread_count);

#pragma omp parallel
    {
        double start_time = wtime();

        int th_id = omp_get_thread_num();
        int remaining_size = (th_id * avg_size) + avg_size;

        std::cout << "Thread id : " << th_id << std::endl;

        if(th_id == thread_count - 1){
            remaining_size = candidates_count[start_vertex];
        }

        //Allocate Memory Buffer

        ui *idx = new ui[max_depth];
        ui *idx_count = new ui[max_depth];
        ui *embedding = new ui[max_depth];
        bool *visited_vertices = new bool[data_graph->getVerticesCount()];
        std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
        VertexID **valid_candidate = new ui *[max_depth];

        ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
        for (ui i = 0; i < max_depth; ++i) {
            valid_candidate[i] = new VertexID[max_candidate_count];
        }

        int cur_depth = 0;

        idx[cur_depth] = 0;
        idx_count[cur_depth] = remaining_size - (th_id * avg_size);
        std::copy(candidates[start_vertex] + (th_id * avg_size), candidates[start_vertex] + remaining_size,
                  valid_candidate[cur_depth]);

        while (true) {
            while (idx[cur_depth] < idx_count[cur_depth]) {
                VertexID u = order[cur_depth];
                VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
                embedding[u] = v;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1) {
                    embedding_cnt_array[th_id][0] += 1;
                    visited_vertices[v] = false;
                    //Enumerate::printMatch(embedding, query_graph->getVerticesCount());
                    if (embedding_cnt_array[th_id][0] >= thread_output_limit_num) {
                        goto EXIT;
                    }
                } else {
                    call_count += 1;
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    Enumerate::generateValidCandidatesWithCandidateCSR(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                            visited_vertices, tree, order, candidates, candidates_count, candidate_offset, candidate_csr);
                }
            }

            cur_depth -= 1;
            if (cur_depth < 0)
                break;
            else
                visited_vertices[embedding[order[cur_depth]]] = false;
        }

        double end_time = wtime();
        thread_wise_time[th_id] = end_time - start_time;

        // Release the buffer.
        EXIT:
        delete[] idx;
        delete[] idx_count;
        delete[] embedding;
        delete[] visited_vertices;
        for (ui i = 0; i < max_depth; ++i) {
            delete[] valid_candidate[i];
        }

        delete[] valid_candidate;

    }

    std::cout << "Thread wise time spent for " << thread_count << " threads" << std::endl;

    for(ui i = 0; i < thread_count; i++){
        std::cout << "Thread ID : " << i << " Spent Time : " << thread_wise_time[i] << std::endl;
    }

    return embedding_cnt_array;
}

void ParallelEnumeration::writeResult(const std::string& file_path, int& thread_count, size_t &call_count, size_t &embedding_count){
    std::ofstream outputfile;
    outputfile.open(file_path, std::ios::app);

    outputfile << "----------------------------------------" << std::endl;
    outputfile << std::endl;

    outputfile << "Thread Count : " << thread_count << std::endl;
    outputfile << "Embedding Count : " << embedding_count << std::endl;
    outputfile << "Operation Count : " << call_count << std::endl;
    outputfile << std::endl;

    outputfile << "----------------------------------------" << std::endl;

    outputfile.flush();
    outputfile.close();
}


ui* ParallelEnumeration::explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
               TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count, int &thread_count){

    std::cout << " ################## explore parallel ##################" << std::endl;

    ui* embedding_cnt_array = new ui[thread_count];

    for(ui i = 0; i < NUM_THREADS; i++){
        embedding_cnt_array[i] = 0;
    }

    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    int avg_size = candidates_count[start_vertex] / thread_count;

    // Allocate the memory buffer.

    omp_set_num_threads(thread_count);

#pragma omp parallel
    {

        int th_id = omp_get_thread_num();
        int remaining_size = (th_id * avg_size) + avg_size;

        std::cout << "Thread id : " << th_id << std::endl;

        if(th_id == NUM_THREADS - 1){
            remaining_size = candidates_count[start_vertex];
        }

        ui *idx = new ui[max_depth];
        ui *idx_count = new ui[max_depth];
        ui *embedding = new ui[max_depth];
        bool *visited_vertices = new bool[data_graph->getVerticesCount()];
        std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
        VertexID **valid_candidate = new ui *[max_depth];

        ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
        for (ui i = 0; i < max_depth; ++i) {
            valid_candidate[i] = new VertexID[max_candidate_count];
        }

        int cur_depth = 0;

        idx[cur_depth] = 0;
        idx_count[cur_depth] = remaining_size - (th_id * avg_size);
        std::copy(candidates[start_vertex] + (th_id * avg_size), candidates[start_vertex] + remaining_size,
                  valid_candidate[cur_depth]);

        while (true) {
            while (idx[cur_depth] < idx_count[cur_depth]) {
                VertexID u = order[cur_depth];
                VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
                embedding[u] = v;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1) {
                    embedding_cnt_array[th_id] += 1;
                    visited_vertices[v] = false;
                    //Enumerate::printMatch(embedding, query_graph->getVerticesCount());
                    if (embedding_cnt_array[th_id] >= thread_output_limit_num) {
                        goto EXIT;
                    }
                } else {
                    call_count += 1;
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    Enumerate::generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                            visited_vertices, tree, order, candidates, candidates_count);
                }
            }

            cur_depth -= 1;
            if (cur_depth < 0)
                break;
            else
                visited_vertices[embedding[order[cur_depth]]] = false;
        }

        // Release the buffer.
        EXIT:
        delete[] idx;
        delete[] idx_count;
        delete[] embedding;
        delete[] visited_vertices;
        for (ui i = 0; i < max_depth; ++i) {
            delete[] valid_candidate[i];
        }

        delete[] valid_candidate;

    }

    return embedding_cnt_array;
}