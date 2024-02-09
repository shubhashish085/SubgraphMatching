#ifndef BACKTRACKING_H
#define BACKTRACKING_H

#include "graph.h"
#include <vector>
#include <utility>

class AlgorithmStore{

public:
    static bool DFSTraversal(const Graph *graph, ui start_node, std::vector<std::pair<ui, ui>>& non_tree_edges, std::vector<ui>& matching_order, std::vector<bool>& visited);
    static bool BFSTraversal(const Graph* graph, ui start_node, std::vector<std::pair<ui, ui>>& non_tree_edges, std::vector<ui>& matching_order, std::vector<bool>& visited);
};

#endif