//
// Created by antu on 3/11/24.
//

#include "utilities.h"
#include "Enumeration.h"

/* Tree list implementation */

void Enumerate::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                            ui *order, ui **candidates, ui *candidates_count) {

    //std::cout << " ################## generateValidCandidates ##################" << std::endl;

    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < tree[u].bn_count_; ++j) {
                VertexID u_nbr = tree[u].bn_[j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }

    //std::cout << " ################## end generateValidCandidates ##################" << std::endl;
}

void Enumerate::generateValidCandidatesForRecursive(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                        ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                        ui *order, ui **candidates, ui *candidates_count) {


    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < tree[u].bn_count_; ++j) {
                VertexID u_nbr = tree[u].bn_[j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}


void Enumerate::printMatch(ui* embedding, ui max_depth){

    std::cout << "------------------ Matched -------------------" << std::endl;

    for(ui i = 0; i < max_depth; i++){
        std::cout << "(" << i << "," << embedding[i] << "),";
    }

    std::cout << std::endl;
}


size_t Enumerate::explore(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order, TreeNode *& tree,
                                          size_t output_limit_num, size_t &call_count) {

    std::cout << " ################## explore ##################" << std::endl;

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    call_count = 0;


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
    std::cout << "Max candidate count : " << max_candidate_count << std::endl;
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }



    std::cout << "Candidate count of Start Vertex : " << candidates_count[start_vertex] << std::endl;
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    std::cout << "Entering Loop " << std::endl;

    while (true) {
        //std::cout << " idx[cur_depth] : " << idx[cur_depth] << " idx_count[cur_depth] : " << idx_count[cur_depth] << std::endl;
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            //std::cout << " v : " << v << std::endl;
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;
            call_count++;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                //printMatch(embedding, query_graph->getVerticesCount());
                if (embedding_cnt >= output_limit_num) {
                    std::cout << "Output Limit Exceeded" << std::endl;
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, tree, order, candidates, candidates_count);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else if (cur_depth == 0){
            std::cout << "Candidate covered : " << idx[cur_depth] << " of candidate count : " << idx_count[cur_depth] << " at level cur_depth : " << cur_depth << std::endl;
            visited_vertices[embedding[order[cur_depth]]] = false;
        }else {
            visited_vertices[embedding[order[cur_depth]]] = false;
        }
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

    std::cout << " Total Embedding Count : " << embedding_cnt << std::endl;
    std::cout << "Total Operation Count : " << call_count << std::endl;

    return embedding_cnt;
}

void Enumerate::exploreRecursiveFashion(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *embedding,
                               ui curr_depth, ui max_depth, ui *order, ui* idx, ui* idx_count, bool* visited_vertices, VertexID **valid_candidate,
                               TreeNode *& tree, size_t &embedding_count, size_t &call_count){

    //std::cout << " Current Depth : " << curr_depth << std::endl;

    if(curr_depth == query_graph ->getVerticesCount() - 1){
        VertexID u = order[curr_depth];
        VertexID v = valid_candidate[curr_depth][idx[curr_depth]];
        embedding[u] = v;
        //printMatch(embedding, max_depth);
        embedding_count++;
        return;

    }else{

        VertexID u = order[curr_depth];
        VertexID v = valid_candidate[curr_depth][idx[curr_depth]];
        std::cout << " u: " << u << " v: " << v << std::endl;
        embedding[u] = v;
        visited_vertices[v] = true;
        idx[curr_depth]++;


        call_count += 1;
        curr_depth += 1;
        idx[curr_depth] = 0;
        generateValidCandidatesForRecursive(data_graph, curr_depth, embedding, idx_count, valid_candidate,
                                                visited_vertices, tree, order, candidates, candidates_count);

        //std::cout << " Valid Candidate Count : " << idx_count[curr_depth] << std::endl;

        for(ui i = 0; (curr_depth <= max_depth - 1) && (i < idx_count[curr_depth]); i++){
            exploreRecursiveFashion(data_graph, query_graph, candidates, candidates_count, embedding,
                                    curr_depth, max_depth, order, idx, idx_count, visited_vertices, valid_candidate,
                                    tree, embedding_count, call_count);
            idx[curr_depth]++;
        }

    }

    curr_depth -= 1;
    if(curr_depth == 0){
        visited_vertices[embedding[order[curr_depth]]] = false;
        return;
    }else {
        visited_vertices[embedding[order[curr_depth]]] = false;
    }

}

void Enumerate::exploreWithoutCandidate(const Graph *data_graph, const Graph *query_graph, ui *order, ui *embedding,
                                        ui curr_depth, bool* visited_vertices, TreeNode *& tree, size_t &embedding_count, size_t &call_count) {

    call_count++;

    if(curr_depth > query_graph->getVerticesCount() - 1){
        embedding_count++;
        //printMatch(embedding, query_graph->getVerticesCount());
        return;
    }

    VertexID u = order[curr_depth];

    //determining v_candidates and candidate_count
    VertexID * v_candidates = new VertexID [data_graph->getGraphMaxLabelFrequency()];
    ui candidate_count = 0;

    if(curr_depth == 0 && tree[u].bn_count_ == 0){

        for(ui i = 0; i < data_graph -> getVerticesCount(); i++){
            if((data_graph ->getVertexLabel(i) == query_graph ->getVertexLabel(u) && data_graph ->getVertexDegree(i) >= query_graph ->getVertexDegree(u))){
                v_candidates[candidate_count++] = i;
            }
        }
    }

    for(ui i = 0; i < tree[u].bn_count_; i++){
        VertexID prev_v = embedding[tree[u].bn_[i]];
        ui nbrs_count = 0;
        VertexID * prev_v_nbrs = data_graph ->getVertexNeighbors(prev_v, nbrs_count);

        ui valid_nbrs_count = 0;
        for(ui j = 0; j < nbrs_count; j++){
            if((data_graph ->getVertexLabel(prev_v_nbrs[j]) == query_graph ->getVertexLabel(u) && data_graph ->getVertexDegree(prev_v_nbrs[j]) >= query_graph ->getVertexDegree(u))){
                prev_v_nbrs[valid_nbrs_count++] = prev_v_nbrs[j];
            }
        }

        //set intersection
        if(i == 0){
            Utilities::set_copy(prev_v_nbrs, valid_nbrs_count, v_candidates, candidate_count);
        }else{
            Utilities::set_intersection(prev_v_nbrs, valid_nbrs_count, v_candidates, candidate_count);
        }

    }

    //std::cout << "Candidates of " << u << " : " ;
    for(ui i = 0; i < candidate_count; i++){
        std::cout << v_candidates[i] << ", ";
    }
    std::cout << std::endl;


    for(ui i = 0; i < candidate_count; i++){

        VertexID v = v_candidates[i];

        if(!visited_vertices[v] && (data_graph ->getVertexLabel(v) == query_graph ->getVertexLabel(u) &&
            data_graph ->getVertexDegree(v) >= query_graph ->getVertexDegree(u))){
            embedding[u] = v;
            visited_vertices[v] = true;
            curr_depth++;
            exploreWithoutCandidate(data_graph, query_graph, order, embedding,
                    curr_depth, visited_vertices, tree, embedding_count, call_count);
            curr_depth--;
            visited_vertices[v] = false;
        }
    }

    if(curr_depth <= 0){
        return;
    }

}



