

#include "graph.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>


void Graph::BuildReverseIndex() {
    reverse_index = new ui[vertices_count];
    reverse_index_offsets= new ui[labels_count + 1];
    reverse_index_offsets[0] = 0;

    ui total = 0;
    for (ui i = 0; i < labels_count; ++i) {
        reverse_index_offsets[i + 1] = total;
        total += labels_frequency[i];
    }

    for (ui i = 0; i < vertices_count; ++i) {
        LabelID label = labels[i];
        reverse_index[reverse_index_offsets[label + 1]++] = i;
    }
}

void Graph::loadGraphFromFileWithEdge(const std::string& file_path){
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    std::string input_line;
    ui label = 0;

    std::cout << "Reading File............ " << std::endl;

    while (std::getline(infile, input_line)) {

        std::cout << " Input Line : " << input_line << std::endl;

        if(input_line.rfind("#", 0) == 0) {
            if(input_line.rfind("# Nodes", 0) == 0) {
                std::stringstream ss(input_line);
                std::string token;
                int count = 0;
                while (!ss.eof()) {
                    std::getline(ss, token, ' ');
                    if (!(token.rfind("#", 0) == 0 || token.rfind("Nodes:", 0) == 0 || token.rfind("Edges:", 0) == 0)) {
                        if (count == 0) {
                            vertices_count = stoi(token);
                            std::cout << " Vertex Count : " << vertices_count << std::endl;
                            degrees = new ui[vertices_count];
                            std::fill(degrees, degrees + vertices_count, 0);
                            count = 1;
                        } else {
                            edges_count = stoi(token);
                            count = 0;
                        }
                        std::cout << "Vertices Count : " << vertices_count << " Edges Count : " << edges_count << std::endl;
                    }
                }
            }
        }else{
            std::stringstream ss(input_line);
            std::string token;
            int count = 0, index;
            while (!ss.eof()) {
                std::getline(ss, token, ' ');
                index = stoi(token);
                degrees[index] += 1;
            }

        }
    }

    infile.close();

    std::ifstream input_file(file_path);

    offsets = new ui[vertices_count +  1];
    offsets[0] = 0;

    neighbors = new VertexID[edges_count * 2];
    labels = new LabelID[vertices_count];
    neighborhood_label_count = new std::unordered_map<LabelID, ui>[vertices_count];
    labels_count = 0;
    max_degree = 0;

    LabelID max_label_id = 0, begin_vtx_label, end_vtx_label;
    std::vector<ui> neighbors_offset(vertices_count, 0);// used for adjust neighbors with offset

    for(ui id = 0; id < vertices_count; id++){
        labels[id] = label;
        offsets[id + 1] = offsets[id] + degrees[id];

        if (degrees[id] > max_degree) {
            max_degree = degrees[id];
        }

        if (labels_frequency.find(label) == labels_frequency.end()) {
            labels_frequency[label] = 0;
            if (label > max_label_id)
                max_label_id = label;
        }

        labels_frequency[label] += 1;
    }

    ui line_count = 0, count = 0;

    while (std::getline(input_file, input_line)) {

        if(line_count < 4){

            line_count++;

        }else{ // Read edge.
            VertexID begin;
            VertexID end;

            std::stringstream ss(input_line);
            std::string token;
            count = 0;
            while (!ss.eof()) {
                std::getline(ss, token, ' ');
                if(count == 0){
                    begin = stoi(token);
                    count = 1;
                }else{
                    end = stoi(token);
                    count = 0;
                }
            }

            ui offset = offsets[begin] + neighbors_offset[begin]; // adjusting the index of neighbor in neighbors array
            neighbors[offset] = end;

            offset = offsets[end] + neighbors_offset[end]; // adjusting the index of neighbor in neighbors array
            neighbors[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;

            if(neighborhood_label_count[begin].find(labels[end]) == neighborhood_label_count[end].end()){
                neighborhood_label_count[begin][labels[end]] = 0;
            }
            neighborhood_label_count[begin][labels[end]] += 1;

            if(neighborhood_label_count[end].find(labels[begin]) == neighborhood_label_count[begin].end()){
                neighborhood_label_count[end][labels[begin]] = 0;
            }
            neighborhood_label_count[end][labels[begin]] += 1;

        }
    }

    std::cout << std::endl;

    input_file.close();
    labels_count = (ui)labels_frequency.size() > (max_label_id + 1) ? (ui)labels_frequency.size() : max_label_id + 1;

    for (auto element : labels_frequency) {
        if (element.second > max_label_frequency) {
            max_label_frequency = element.second;
        }
    }

    for (ui i = 0; i < vertices_count; ++i) {
        std::sort(neighbors + offsets[i], neighbors + offsets[i + 1]); // sorting the neighbors of every vertex
    }

    BuildReverseIndex();

    printGraphData();
}

void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ." << std::endl;
        exit(-1);
    }

    char type;
    infile >> type >> vertices_count >> edges_count;
    offsets = new ui[vertices_count +  1];
    offsets[0] = 0;

    neighbors = new VertexID[edges_count * 2];
    labels = new LabelID[vertices_count];
    neighborhood_label_count = new std::unordered_map<LabelID, ui>[vertices_count];
    labels_count = 0;
    max_degree = 0;

    LabelID max_label_id = 0, begin_vtx_label, end_vtx_label;
    std::vector<ui> neighbors_offset(vertices_count, 0);// used for adjust neighbors with offset

    while (infile >> type) {
        if (type == 'v') { // Reading vertex.
            VertexID id;
            LabelID  label;
            ui degree;
            infile >> id >> label >> degree;

            labels[id] = label;
            offsets[id + 1] = offsets[id] + degree;

            if (degree > max_degree) {
                max_degree = degree;
            }

            if (labels_frequency.find(label) == labels_frequency.end()) {
                labels_frequency[label] = 0;
                if (label > max_label_id)
                    max_label_id = label;
            }

            labels_frequency[label] += 1;
        }
        else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            infile >> begin >> end;

            ui offset = offsets[begin] + neighbors_offset[begin]; // adjusting the index of neighbor in neighbors array
            neighbors[offset] = end;

            offset = offsets[end] + neighbors_offset[end]; // adjusting the index of neighbor in neighbors array
            neighbors[offset] = begin;

            neighbors_offset[begin] += 1;
            neighbors_offset[end] += 1;

            if(neighborhood_label_count[begin].find(labels[end]) == neighborhood_label_count[end].end()){
                neighborhood_label_count[begin][labels[end]] = 0;
            }
            neighborhood_label_count[begin][labels[end]] += 1;

            if(neighborhood_label_count[end].find(labels[begin]) == neighborhood_label_count[begin].end()){
                neighborhood_label_count[end][labels[begin]] = 0;
            }
            neighborhood_label_count[end][labels[begin]] += 1;

        }
    }

    /*std::cout << " --------- Neighbors " << std::endl;

    for(ui i = 0; i < edges_count * 2; i++){
        std::cout << " " << neighbors[i] ;
    }

    std::cout << std::endl;

    std::cout << " --------- Offsets " << std::endl;

    for(ui i = 0; i < vertices_count + 1; i++){
        std::cout << " " << offsets[i] ;
    }*/

    std::cout << std::endl;

    infile.close();
    labels_count = (ui)labels_frequency.size() > (max_label_id + 1) ? (ui)labels_frequency.size() : max_label_id + 1;

    for (auto element : labels_frequency) {
        if (element.second > max_label_frequency) {
            max_label_frequency = element.second;
        }
    }

    for (ui i = 0; i < vertices_count; ++i) {
        std::sort(neighbors + offsets[i], neighbors + offsets[i + 1]); // sorting the neighbors of every vertex
    }

    BuildReverseIndex();

}

void Graph::setMatchingOrderIndex(std::vector<ui> matching_order){

    matching_order_idx = new ui[vertices_count];

    for(ui i = 0; i < vertices_count; i++){
        matching_order_idx[matching_order[i]] = i;
    }

}

void Graph::printGraphMetaData() {
    std::cout << "|V|: " << vertices_count << ", |E|: " << edges_count << ", |\u03A3|: " << labels_count << std::endl;
    std::cout << "Max Degree: " << max_degree << ", Max Label Frequency: " << max_label_frequency << std::endl;
}