//
// Created by kars1 on 5/23/24.
//

#ifndef SUBGRAPHMATCHING_BUILDTABLE_H
#define SUBGRAPHMATCHING_BUILDTABLE_H

#include "graph.h"

class BuildTable {

public:
    static size_t computeMemoryCostInBytes(const Graph *query_graph, const Graph *data_graph, ui *candidates_count);

};


#endif //SUBGRAPHMATCHING_BUILDTABLE_H
