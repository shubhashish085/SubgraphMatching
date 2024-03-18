//
// Created by antu on 3/15/24.
//

#include "utilities.h"

void Utilities::set_copy(VertexID *l_array, ui &l_count, VertexID *valid_candidate, ui &candidate_count) {

    candidate_count = 0;

    for(ui i = 0; i < l_count; i++){
        valid_candidate[candidate_count++] = l_array[i];
    }

    return;
}

void Utilities::set_intersection(VertexID *r_array, ui &r_count,
                                 VertexID *valid_candidate, ui &candidate_count) {

    if(candidate_count == 0){
        return;
    }
    VertexID * l_array = new VertexID[candidate_count];
    ui  l_count = candidate_count;

    for(ui i = 0; i < l_count; i++){
        l_array[i] = valid_candidate[i];
    }

    candidate_count = 0;

    for(ui i = 0; i < l_count; i++){
        bool shouldBeInserted = false;
        for(ui j = 0; j < r_count; j++){
            if(l_array[i] == r_array[j]){
                shouldBeInserted = true;
                break;
            }
        }

        if(shouldBeInserted){
            valid_candidate[candidate_count++] = r_array[i];
        }
    }

    delete[] l_array;

    return;
}
