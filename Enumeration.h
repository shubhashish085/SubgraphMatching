//
// Created by antu on 3/11/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_ENUMERATE_H
#define SUBGRAPHMATCHINGMAIN_ENUMERATE_H

#include "graph.h"
#include "types.h"

class Enumerate {

public:
    static size_t explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                          TreeNode *& tree, size_t output_limit_num, size_t &call_count);

    static void generateValidCandidates(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count);


};


#endif //SUBGRAPHMATCHINGMAIN_ENUMERATE_H
