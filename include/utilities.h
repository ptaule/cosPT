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



inline short int sum_vectors(
        const short int labels[],
        size_t n_vecs,
        const short int sum_table[][N_CONFIGS]
        )
{
    if (n_vecs == 1) return labels[0];

    short int result = labels[0];

    for (size_t i = 1; i < n_vecs; ++i) {
        if (labels[i] == ZERO_LABEL) continue;
        result = sum_table[result][labels[i]];
    }

    // If DEBUG==true, check that sum is an appropriate vector configuration,
    // i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is
    // an element of (0,1)
#if DEBUG >= 1
    short int res_coeffs[N_COEFFS];
    label2config(result,res_coeffs,N_COEFFS);
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        short int c = res_coeffs[i];
        if (!(c == -1 || c == 0 || c == 1))
            warning("Sum of vectors does not correspond to an appropriate "
                    "configuration.");
    }
    short int c = res_coeffs[N_COEFFS - 1];
    if (!(c == 0 || c == 1))
        warning("Sum of vectors does not correspond to an appropriate "
                "configuration.");
#endif

    return result;
}



void compute_sum_table(short int sum_table[][N_CONFIGS]);

short int zero_label();

bool is_fundamental(short int label);
bool unique_elements(const short int array[],size_t length, short int skip);

void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width);


#endif /* ifndef UTILITIES_H */
