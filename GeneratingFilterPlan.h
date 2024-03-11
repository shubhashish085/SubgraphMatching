//
// Created by antu on 3/10/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_GENERATINGFILTERPLAN_H
#define SUBGRAPHMATCHINGMAIN_GENERATINGFILTERPLAN_H

#include "graph.h"
#include "types.h"

class GeneratingFilterPlan {

public:
    static void generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                      VertexID *&order, int &level_count, ui *&level_offset);

private:
    static VertexID selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph);
};


#endif //SUBGRAPHMATCHINGMAIN_GENERATINGFILTERPLAN_H
