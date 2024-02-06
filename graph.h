
#ifndef GRAPH_H
#define GRAPH_H

#include "types.h"

class Graph{

private:

    ui vertices_count;
    ui edges_count;
    ui labels_count;
    ui max_degree;
    ui max_label_frequency;

    ui* offsets;
    VertexID * neighbors;
    LabelID* labels;

    std::unordered_map<LabelID, ui> labels_frequency;

public:

    Graph() {

        vertices_count = 0;
        edges_count = 0;
        labels_count = 0;
        max_degree = 0;
        max_label_frequency = 0;

        offsets = NULL;
        neighbors = NULL;
        labels = NULL;

        labels_frequency.clear();

    }

    ~Graph() {
        delete[] offsets;
        delete[] neighbors;
        delete[] labels;
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void printGraphMetaData();

public:
    const ui getLabelsCount() const {
        return labels_count;
    }

    const ui getVerticesCount() const {
        return vertices_count;
    }

    const ui getEdgesCount() const {
        return edges_count;
    }

    const ui getGraphMaxDegree() const {
        return max_degree;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency;
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency.find(label) == labels_frequency.end() ? 0 : labels_frequency.at(label);
    }

    const LabelID getVertexLabel(const VertexID id) const {
        return labels[id];
    }

    const ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id]; // used for neighbor count
        return neighbors + offsets[id];
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

};


#endif