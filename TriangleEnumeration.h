//
// Created by antu on 6/14/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H
#define SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H

#include "graph.h"
#include "types.h"

class TriangleEnumeration {

public:
    static size_t enumerateTriangles(const Graph *data_graph, size_t output_limit_num, size_t &call_count);
    static void buildGraph();



};


#endif //SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H
