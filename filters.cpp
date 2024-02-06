
#include "filters.h"
#include <vector>
#include <algorithm>


bool
Filters::LDFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        LabelID label = query_graph->getVertexLabel(i);
        ui degree = query_graph->getVertexDegree(i);

        ui data_vertex_num;
        const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);

        for (ui j = 0; j < data_vertex_num; ++j) {
            ui data_vertex = data_vertices[j];
            if (data_graph->getVertexDegree(data_vertex) >= degree) {
                candidates[i][candidates_count[i]++] = data_vertex;
            }
        }

        if (candidates_count[i] == 0) {
            return false;
        }
    }

    return true;
}

bool
Filters::NLFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        computeCandidateWithNLF(data_graph, query_graph, query_vertex, candidates_count[query_vertex], candidates[query_vertex]);

        if (candidates_count[query_vertex] == 0) {
            return false;
        }
    }

    return true;
}




void Filters::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    printf("############## Allocating Buffer ##############");
    ui query_vertex_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency(); // candidate number of any vertex can be at most the maximum label frequency of
    //  any label in the data graph

    candidates_count = new ui[query_vertex_num];// storing the candidate count for every vertex in the query graph
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui*[query_vertex_num]; // storing the candidates for every vertex in the query graph

    for (ui i = 0; i < query_vertex_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}