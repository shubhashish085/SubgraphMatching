
#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
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
    ui* matching_order_idx;
    LabelID* labels;

    std::unordered_map<LabelID, ui> * neighborhood_label_count;
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
        neighborhood_label_count = NULL;

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
    void setMatchingOrderIndex(std::vector<ui> matching_order);

public:

    const ui* getOffsets() const {
        return offsets;
    }

    const ui getLabelsCount() const {
        return labels_count;
    }

    const ui getVerticesCount() const {
        return vertices_count;
    }

    const ui getEdgesCount() const {
        return edges_count;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets[id + 1] - offsets[id];
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

    ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id]; // used for neighbor count
        return neighbors + offsets[id];
    }

    ui * getMatchingOrderIndex() const {
        return matching_order_idx;
    }

    std::unordered_map<LabelID, ui> * getNeighborhoodLabelCount(){
        return neighborhood_label_count;
    }

    const std::vector<VertexID> getVerticesOfLabel(LabelID labelId) const {

        std::vector<VertexID> candidate_vertices;

        for(VertexID i = 0; i < vertices_count; i++){
            if(labels[i] == labelId){
                candidate_vertices.push_back(i);
            }
        }

        return candidate_vertices;
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