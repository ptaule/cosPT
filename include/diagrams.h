/*
   diagrams.h

   Created by Petter Taule on 27.03.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef DIAGRAMS_H
#define DIAGRAMS_H

#include "constants.h"
#include "tables.h"

// Struct for labeling different diagrams
typedef struct {
    short int    m;
    short int    l;
    short int    r;
    short int    diagram_factor;     /* Topological multiplicative diagram factor      */
    short int    n_rearrangements;   /* Number of rearrangements of loop momenta       */
    short int    n_sign_configs;     /* Number of sign flips                           */
    short int**  rearrangements;     /* Table of rearrangements                        */
    bool**       sign_configs;       /* Table of sign flips, true <-> +1, false <-> -1 */
    short int**  kernel_indices_l;
    short int**  kernel_indices_r;
    short int*** argument_configs_l;
    short int*** argument_configs_r;
} diagram_t;

void initialize_diagrams(diagram_t diagrams[]);
void diagrams_gc(diagram_t diagrams[]);

void print_diagram(const diagram_t* diagram);

#endif /* ifndef DIAGRAMS_H */
