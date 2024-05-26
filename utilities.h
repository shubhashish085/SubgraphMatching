//
// Created by antu on 3/15/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_UTILITIES_H
#define SUBGRAPHMATCHINGMAIN_UTILITIES_H

#include "types.h"

class Utilities {

public:
    static void set_copy (VertexID* l_array, ui &l_count, VertexID* valid_candidate, ui &candidate_count);

    static void set_intersection_tp(VertexID * result, ui l_length, VertexID * r_array, ui r_length, ui& set_ints_length);

    static void set_intersection_with_maximum_bound(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
            ui &set_ints_length, ui max);

    static void set_intersection (VertexID* r_array, ui &r_count, VertexID* valid_candidate, ui &candidate_count);

    static int binary_search (ui* array, ui low_idx, ui high_idx, ui element);

};


#endif //SUBGRAPHMATCHINGMAIN_UTILITIES_H
