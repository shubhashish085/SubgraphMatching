//
// Created by antu on 3/18/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_PARALLELENUMERATION_H
#define SUBGRAPHMATCHINGMAIN_PARALLELENUMERATION_H


#include "types.h"
#include "graph.h"


class ParallelEnumeration {

public:
    static void writeResult(const std::string& file_path, int& thread_count, size_t &call_count, size_t &embedding_count);
    static ui* explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                            TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count, int &thread_count);
    static ui** exploreWithPadding(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                 TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count, int &thread_count);

};


#endif //SUBGRAPHMATCHINGMAIN_PARALLELENUMERATION_H
