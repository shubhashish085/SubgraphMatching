//
// Created by antu on 6/14/24.
//

#include <map>
#include <vector>
#include "backtracking.h"
#include "TriangleEnumeration.h"
#include "utilities.h"

size_t TriangleEnumerate::enumerateTriangles(const Graph* data_graph, std::vector<std::pair<VertexID, VertexID>>& edges, ui* result_array,
                                               size_t &output_limit_num, size_t &call_count) {

    size_t  triangle_count = 0;

    //BFS at all node
    std::map<VertexID, std::vector<VertexID>> u_w_map;
    std::map<VertexID, std::vector<VertexID>>::iterator u_w_map_iter;
    ui u_level, v_level, distance = 2, edge_sg_count = 0, edge_dg_count = 0, offset;
    VertexID small_idx_vtx, large_idx_vtx;


    ui vertex_num = data_graph->getVerticesCount();
    ui* two_hop_path_count_sg = new ui[vertex_num];
    ui* two_hop_path_count_dg = new ui[vertex_num];
    ui* distance_array = new ui[vertex_num];
    VertexID* explored_vertices = new VertexID[vertex_num];
    bool* visited = new bool[vertex_num];
    int* level = new int[vertex_num];

    ui* neighbor_offsets_sg = new ui [vertex_num];
    ui* neighbor_offsets_dg = new ui [vertex_num];
    ui* offsets_sg = new ui[vertex_num + 1];
    ui* offsets_dg = new ui[vertex_num + 1];
    ui* degree_sg = new ui[vertex_num];
    ui* degree_dg = new ui[vertex_num];

    for(ui i = 0; i < vertex_num; i++){
        level[i] = -1;
        result_array[i] = 0;
        two_hop_path_count_sg[i] = 0;
        two_hop_path_count_dg[i] = 0;
        distance_array[i] = 0;
        degree_sg[i] = 0;
        degree_dg[i] = 0;
        offsets_sg[i] = 0;
        offsets_dg[i] = 0;
        neighbor_offsets_sg[i] = 0;
        neighbor_offsets_dg[i] = 0;
    }

    offsets_sg[vertex_num] = 0;
    offsets_dg[vertex_num] = 0;

    for(ui i = 0; i < vertex_num; i++){
        if(level[i] == -1){
            AlgorithmStore::bfsTraversal(data_graph, i, level);
        }
    }


    std::cout << std::endl;

    std::vector<std::pair<VertexID, VertexID>>::iterator iter = edges.begin();

    //std::cout << "2. Started Edge Separation " << std::endl;

    //level wise edge separation and graph building
    for (iter; iter != edges.end(); iter++){

        u_level = level[iter->first];
        v_level = level[iter->second];

        small_idx_vtx = iter->first < iter -> second ? iter -> first: iter -> second;
        large_idx_vtx = iter->first < iter -> second ? iter -> second: iter -> first;

        if (u_level == v_level) {

            //std::cout << "Same level vertices : " << small_idx_vtx << ", " << large_idx_vtx << std::endl;

            u_w_map_iter = u_w_map.find(small_idx_vtx);
            if(u_w_map_iter != u_w_map.end()){
                (u_w_map_iter->second).push_back(large_idx_vtx);
            } else{
                std::vector<VertexID> vtr;
                vtr.push_back(large_idx_vtx);
                u_w_map[small_idx_vtx] = vtr;
            }
            degree_sg[iter->first] += 1;
            degree_sg[iter->second] += 1;
            edge_sg_count++;


        } else {
            //std::cout << "Different level vertices : " << small_idx_vtx << ", " << large_idx_vtx << std::endl;
            degree_dg[iter->first] += 1;
            degree_dg[iter->second] += 1;
            edge_dg_count++;
        }

    }


    for(ui i = 0; i < vertex_num; i++){
        offsets_sg[i + 1] = offsets_sg[i] + degree_sg[i];
        offsets_dg[i + 1] = offsets_dg[i] + degree_dg[i];
    }

    VertexID* neighbors_sg = new VertexID [edge_sg_count * 2];
    VertexID* neighbors_dg = new VertexID [edge_dg_count * 2];

    //std::cout << "Setting Neighbors" << std::endl;
    for (iter = edges.begin(); iter != edges.end(); iter++){

        u_level = level[iter->first];
        v_level = level[iter->second];

        if(u_level == v_level){
            offset = offsets_sg[iter->first] + neighbor_offsets_sg[iter->first];
            neighbors_sg[offset] = iter->second;

            offset = offsets_sg[iter->second] + neighbor_offsets_sg[iter->second];
            neighbors_sg[offset] = iter->first;

            neighbor_offsets_sg[iter->first] += 1;
            neighbor_offsets_sg[iter->second] += 1;
        }else {
            offset = offsets_dg[iter->first] + neighbor_offsets_dg[iter->first];
            neighbors_dg[offset] = iter->second;

            offset = offsets_dg[iter->second] + neighbor_offsets_dg[iter->second];
            neighbors_dg[offset] = iter->first;

            neighbor_offsets_dg[iter->first] += 1;
            neighbor_offsets_dg[iter->second] += 1;
        }
    }


    Graph* graph_sg = new Graph(vertex_num, edge_sg_count, offsets_sg, neighbors_sg);
    Graph* graph_dg = new Graph(vertex_num, edge_dg_count, offsets_dg, neighbors_dg);

    for (auto map_iter = u_w_map.begin(); map_iter != u_w_map.end(); map_iter++){

        AlgorithmStore::bfsTraversalWithDistance(graph_sg, map_iter->first, distance_array, explored_vertices, visited, two_hop_path_count_sg, distance);

        for (auto vtr_iter = (map_iter->second).begin(); vtr_iter != (map_iter->second).end(); vtr_iter++){
            triangle_count += two_hop_path_count_sg[*vtr_iter];
            two_hop_path_count_sg[*vtr_iter] = 0;
        }

    }

    triangle_count /= 3;


    for (auto map_iter = u_w_map.begin(); map_iter != u_w_map.end(); map_iter++){

        AlgorithmStore::bfsTraversalWithDistance(graph_dg, map_iter->first, distance_array, explored_vertices, visited, two_hop_path_count_dg, distance);

        for (auto vtr_iter = (map_iter->second).begin(); vtr_iter != (map_iter->second).end(); vtr_iter++){

            triangle_count += two_hop_path_count_dg[*vtr_iter];
            two_hop_path_count_dg[*vtr_iter] = 0;
        }

    }

    return triangle_count;
}

size_t TriangleEnumerate::enumerateTrianglesBySetIntersection(const Graph* data_graph, std::vector<std::pair<VertexID, VertexID>>& edges, ui* result_array,
                                             size_t &output_limit_num, size_t &call_count) {

    size_t  triangle_count = 0;

    //BFS at all node
    std::map<VertexID, std::vector<VertexID>> u_w_map;
    std::map<VertexID, std::vector<VertexID>>::iterator u_w_map_iter;
    ui u_level, v_level, distance = 2, edge_sg_count = 0, edge_dg_count = 0, offset;
    VertexID small_idx_vtx, large_idx_vtx;


    ui vertex_num = data_graph->getVerticesCount();
    ui* two_hop_path_count_sg = new ui[vertex_num];
    ui* two_hop_path_count_dg = new ui[vertex_num];
    ui* distance_array = new ui[vertex_num];
    VertexID* explored_vertices = new VertexID[vertex_num];
    bool* visited = new bool[vertex_num];
    int* level = new int[vertex_num];

    ui* neighbor_offsets_sg = new ui [vertex_num];
    ui* neighbor_offsets_dg = new ui [vertex_num];
    ui* offsets_sg = new ui[vertex_num + 1];
    ui* offsets_dg = new ui[vertex_num + 1];
    ui* degree_sg = new ui[vertex_num];
    ui* degree_dg = new ui[vertex_num];

    for(ui i = 0; i < vertex_num; i++){
        level[i] = -1;
        result_array[i] = 0;
        two_hop_path_count_sg[i] = 0;
        two_hop_path_count_dg[i] = 0;
        distance_array[i] = 0;
        degree_sg[i] = 0;
        degree_dg[i] = 0;
        offsets_sg[i] = 0;
        offsets_dg[i] = 0;
        neighbor_offsets_sg[i] = 0;
        neighbor_offsets_dg[i] = 0;
    }

    offsets_sg[vertex_num] = 0;
    offsets_dg[vertex_num] = 0;

    for(ui i = 0; i < vertex_num; i++){
        if(level[i] == -1){
            AlgorithmStore::bfsTraversal(data_graph, i, level);
        }
    }

    std::cout << std::endl;

    std::vector<std::pair<VertexID, VertexID>>::iterator iter = edges.begin();

    //std::cout << "2. Started Edge Separation " << std::endl;

    //level wise edge separation and graph building
    for (iter; iter != edges.end(); iter++){

        u_level = level[iter->first];
        v_level = level[iter->second];

        small_idx_vtx = iter->first < iter -> second ? iter -> first: iter -> second;
        large_idx_vtx = iter->first < iter -> second ? iter -> second: iter -> first;

        if (u_level == v_level) {

            //std::cout << "Same level vertices : " << small_idx_vtx << ", " << large_idx_vtx << std::endl;

            u_w_map_iter = u_w_map.find(small_idx_vtx);
            if(u_w_map_iter != u_w_map.end()){
                (u_w_map_iter->second).push_back(large_idx_vtx);
            } else{
                std::vector<VertexID> vtr;
                vtr.push_back(large_idx_vtx);
                u_w_map[small_idx_vtx] = vtr;
            }
            degree_sg[iter->first] += 1;
            degree_sg[iter->second] += 1;
            edge_sg_count++;


        } else {
            //std::cout << "Different level vertices : " << small_idx_vtx << ", " << large_idx_vtx << std::endl;
            degree_dg[iter->first] += 1;
            degree_dg[iter->second] += 1;
            edge_dg_count++;
        }

    }

    for(ui i = 0; i < vertex_num; i++){
        offsets_sg[i + 1] = offsets_sg[i] + degree_sg[i];
        offsets_dg[i + 1] = offsets_dg[i] + degree_dg[i];
    }

    VertexID* neighbors_sg = new VertexID [edge_sg_count * 2];
    VertexID* neighbors_dg = new VertexID [edge_dg_count * 2];

    //std::cout << "Setting Neighbors" << std::endl;
    for (iter = edges.begin(); iter != edges.end(); iter++){

        u_level = level[iter->first];
        v_level = level[iter->second];

        if(u_level == v_level){
            offset = offsets_sg[iter->first] + neighbor_offsets_sg[iter->first];
            neighbors_sg[offset] = iter->second;

            offset = offsets_sg[iter->second] + neighbor_offsets_sg[iter->second];
            neighbors_sg[offset] = iter->first;

            neighbor_offsets_sg[iter->first] += 1;
            neighbor_offsets_sg[iter->second] += 1;
        }else {
            offset = offsets_dg[iter->first] + neighbor_offsets_dg[iter->first];
            neighbors_dg[offset] = iter->second;

            offset = offsets_dg[iter->second] + neighbor_offsets_dg[iter->second];
            neighbors_dg[offset] = iter->first;

            neighbor_offsets_dg[iter->first] += 1;
            neighbor_offsets_dg[iter->second] += 1;
        }
    }



    Graph* graph_sg = new Graph(vertex_num, edge_sg_count, offsets_sg, neighbors_sg);
    Graph* graph_dg = new Graph(vertex_num, edge_dg_count, offsets_dg, neighbors_dg);

    ui u_nbr_cnt = 0, v_nbr_cnt = 0, intersection_length;
    VertexID* u_nbrs;
    VertexID* v_nbrs;

    //same level vertices graph
    for (auto map_iter = u_w_map.begin(); map_iter != u_w_map.end(); map_iter++){

        u_nbrs = graph_sg->getVertexNeighbors(map_iter->first, u_nbr_cnt);

        for (auto vtr_iter = (map_iter->second).begin(); vtr_iter != (map_iter->second).end(); vtr_iter++){
            v_nbrs = graph_sg->getVertexNeighbors(*vtr_iter, v_nbr_cnt);
            Utilities::set_intersection_triangle_count(u_nbrs, u_nbr_cnt, v_nbrs, v_nbr_cnt, intersection_length);
            triangle_count += intersection_length;
        }
    }

    triangle_count /= 3;

    //different level vertices graph
    for (auto map_iter = u_w_map.begin(); map_iter != u_w_map.end(); map_iter++){

        u_nbrs = graph_dg->getVertexNeighbors(map_iter->first, u_nbr_cnt);

        for (auto vtr_iter = (map_iter->second).begin(); vtr_iter != (map_iter->second).end(); vtr_iter++){
            v_nbrs = graph_dg->getVertexNeighbors(*vtr_iter, v_nbr_cnt);
            Utilities::set_intersection_triangle_count(u_nbrs, u_nbr_cnt, v_nbrs, v_nbr_cnt, intersection_length);
            triangle_count += intersection_length;
        }
    }

    return triangle_count;
}

