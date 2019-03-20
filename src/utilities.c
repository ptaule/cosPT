/*
   utilities.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <gsl/gsl_matrix.h>

#include "../include/constants.h"
#include "../include/utilities.h"


extern void label2config(short int label, short int config[], size_t size);
extern short int config2label(const short int config[], size_t size);



short int zero_label() {
    const short int coeffs[N_COEFFS] = {};
    return config2label(coeffs, N_COEFFS);
}



bool is_fundamental(short int label){
    // Vector is not fundamental if k is present; this is the case if
    // label >= N_CONFIGS/2
    if (label >= N_CONFIGS/2) return false;

    short int coeffs[N_COEFFS] = {};
    label2config(label,coeffs,N_COEFFS);

    short int num_vecs_present = 0;

    // The last coefficient is for k, hence we can skip this (j < N_COEFFS - 1)
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        if (coeffs[i] != 0) num_vecs_present++;
    }
    return (num_vecs_present == 1);

}



// Are there duplicate elements of the array? Yes, return true;
// no, return false. Do not condiser elements equal to skip.
bool unique_elements(const short int array[], size_t length, short int skip) {
    for (int i = 0; i < length; ++i) {
        short int val = array[i];
        if (val == skip) continue;
        for (int j = i + 1; j < length; ++j) {
            if (array[j] == val) return false;
        }
    }
    return true;
}



short int sum_vectors(const short int labels[], size_t n_vecs) {
    short int temp_coeffs[N_COEFFS] = {};
    short int res_coeffs[N_COEFFS]  = {};

    for (size_t i = 0; i < n_vecs; ++i) {
        if (labels[i] == ZERO_LABEL) continue;
        label2config(labels[i],temp_coeffs,N_COEFFS);
        for (int j = 0; j < N_COEFFS; ++j) {
            res_coeffs[j] += temp_coeffs[j];
        }
    }

    // If DEBUG==true, check that sum is an appropriate vector configuration,
    // i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is
    // an element of (0,1)
#if DEBUG >= 1
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

    return config2label(res_coeffs,N_COEFFS);
}



void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            printf("%g\t",gsl_matrix_get(m,i,j));
        }
        printf("\n");
    }
}
