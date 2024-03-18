//
// Created by antu on 3/15/24.
//

#ifndef SUBGRAPHMATCHINGMAIN_UTILITIES_H
#define SUBGRAPHMATCHINGMAIN_UTILITIES_H

#include "types.h"

class Utilities {

public:
    static void set_copy (VertexID* l_array, ui &l_count, VertexID* valid_candidate, ui &candidate_count);

    static void set_intersection (VertexID* r_array, ui &r_count, VertexID* valid_candidate, ui &candidate_count);

};


#endif //SUBGRAPHMATCHINGMAIN_UTILITIES_H
