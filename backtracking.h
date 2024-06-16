#ifndef BACKTRACKING_H
#define BACKTRACKING_H

#include "graph.h"
#include <vector>
#include <utility>

class AlgorithmStore{

public:
    static bool DFSTraversal(const Graph *graph, ui start_node, std::vector<std::pair<ui, ui>>& non_tree_edges, std::vector<ui>& matching_order, std::vector<bool>& visited);
    static bool BFSTraversal(const Graph* graph, ui start_node, std::vector<std::pair<ui, ui>>& non_tree_edges, std::vector<ui>& matching_order, std::vector<bool>& visited,
                             int* parent_vtr);
    static void bfsTraversal(const Graph *graph, VertexID root_vertex, TreeNode *&tree, VertexID *&bfs_order);
    static void bfsTraversal(const Graph *graph, VertexID root_vertex, ui* level);
    static void bfsTraversalWithDistance(const Graph *graph, VertexID root_vertex, ui* distance_array, ui* path_count, ui distance);
};

#endif