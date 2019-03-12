/*
   utilities.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <gsl/gsl_matrix.h>

#include "../include/constants.h"
#include "../include/utilities.h"



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



short int sum_two_vectors(short int label_a, short int label_b) {
    short int a_coeffs[N_COEFFS]   = {};
    short int b_coeffs[N_COEFFS]   = {};
    short int res_coeffs[N_COEFFS] = {};

    label2config(label_a,a_coeffs,N_COEFFS);
    label2config(label_b,b_coeffs,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        res_coeffs[i] = a_coeffs[i] + b_coeffs[i];
    }

    // If DEBUG==true, check that sum is an appropriate vector configuration, i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is an element of (0,1)
#if DEBUG
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        short int c = res_coeffs[i];
        if (!(c == -1 || c == 0 || c == 1))
            warning_verbose("Sum of vectors with labels (%d,%d) does not "
                    "correspond to an appropriate configuration.", label_a,label_b);
    }
    short int c = res_coeffs[N_COEFFS - 1];
    if (!(c == 0 || c == 1))
        warning_verbose("Sum of vectors with labels (%d,%d) does not "
                "correspond to an appropriate configuration.", label_a,label_b);
#endif

    return config2label(res_coeffs,N_COEFFS);
}



short int sum_vectors(const short int labels[], size_t size) {
    if (size == 1) return labels[0];

    short int result = sum_two_vectors(labels[0],labels[1]);

    for (size_t i = 2; i < size; ++i) {
        if (labels[i] == ZERO_LABEL) continue;
        result = sum_two_vectors(result,labels[i]);
    }
    return result;
}



void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            printf("%g\t",gsl_matrix_get(m,i,j));
        }
        printf("\n");
    }
}
