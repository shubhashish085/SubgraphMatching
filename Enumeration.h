//
// Created by antu on 3/11/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_ENUMERATE_H
#define SUBGRAPHMATCHINGMAIN_ENUMERATE_H

#include "graph.h"
#include "types.h"

class Enumerate {

public:
    static void printMatch(ui* embedding, ui max_depth);

    static void exploreWithoutCandidate(const Graph *data_graph, const Graph *query_graph, ui *order, ui *embedding,
                            ui curr_depth, bool* visited_vertices, TreeNode *& tree);

    static void exploreRecursiveFashion(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *embedding,
                            ui curr_depth, ui max_depth, ui *order, ui* idx, ui* idx_count, bool* visited_vertices, VertexID **valid_candidate,
                            TreeNode *& tree, size_t &embedding_count, size_t &call_count);

    static size_t explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                          TreeNode *& tree, size_t output_limit_num, size_t &call_count);

    static void generateValidCandidates(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count);

    static void generateValidCandidatesForRecursive(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                        ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                        ui *order, ui **candidates, ui *candidates_count);


};


#endif //SUBGRAPHMATCHINGMAIN_ENUMERATE_H
