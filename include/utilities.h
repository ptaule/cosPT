/*
   utilities.h

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdbool.h>
#include <math.h>

#include <gsl/gsl_matrix.h>

#include "constants.h"

// Utility functions

inline void label2config(
        short int label,    // in, label element in [0, N_CONFIG - 1]
        short int config[], // out, array of base 3 digits
        size_t size   // in, size of coefficients
        )
{
    // Q-coefficients are element in {-1,0,1}
    for (size_t i = 0; i < size - 1; ++i) {
        config[i] = (label % 3) - 1;
        label /= 3;
    }

    // k-coefficient (last coefficient) is element in {0,1}

    config[size - 1] = label % 3;
}



inline short int config2label(
        const short int config[], // in, array of base 3 digits
        size_t size         // in, size of coefficients
        )
{
    short int label = 0;

    // For Q-coefficients, add 1 to obtain range {0,1,2}
    for (size_t i = 0; i < size - 1; ++i) {
        label += (config[i] + 1) * pow(3,i);
    }

    // For k-coefficient, do _not_ add 1, range is {0,1}
    label += config[size-1] * pow(3,size - 1);

    return label;
}

void initialize_timesteps(double eta[], double eta_i, double eta_f);

void print_label(short int label);
void print_labels(const short int labels[]);
void print_gsl_matrix(const matrix_t* m, size_t height, size_t width);

short int zero_label();
bool is_fundamental(short int label);
bool unique_elements(const short int array[], short int  length, short int skip);

#endif /* ifndef UTILITIES_H */
