
#ifndef SUBGRAPHMATCHINGMAIN_PRUNINGCONSTRAINTS_H
#define SUBGRAPHMATCHINGMAIN_PRUNINGCONSTRAINTS_H



#include "types.h"
#include "graph.h"
#include <vector>
#include <map>



class PruningConstraints {

public:
    static void buildCandidateMap(const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order, std::map<VertexID, std::vector<VertexID>>& candidate_map);
    static bool isCandidate(VertexID query_vertex, VertexID vertexToBeChecked, std::map<VertexID, std::vector<VertexID>>& candidate_map);
    static void checkingCycle(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                        std::vector<std::tuple<VertexID, VertexID>>& edges, std::map<VertexID, std::vector<VertexID>>& candidate_map);

};

#endif //SUBGRAPHMATCHINGMAIN_PRUNINGCONSTRAINTS_H

