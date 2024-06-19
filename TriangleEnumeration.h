//
// Created by antu on 6/14/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H
#define SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H

#include "graph.h"
#include "types.h"

class TriangleEnumerate {

public:
    static size_t enumerateTriangles(const Graph* data_graph, std::vector<std::pair<VertexID, VertexID>>& edges, ui* result_array,
                                     size_t &output_limit_num, size_t &call_count);
    static size_t enumerateTrianglesBySetIntersection(const Graph* data_graph, std::vector<std::pair<VertexID, VertexID>>& edges, ui* result_array,
                                               size_t &output_limit_num, size_t &call_count);



};


#endif //SUBGRAPHMATCHINGMAIN_TRIANGLEENUMERATION_H
