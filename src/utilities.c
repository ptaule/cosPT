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
extern short int sum_vectors(const short int labels[], size_t n_vecs, const short int sum_table[][N_CONFIGS]);



static inline short int sum_two_vectors(short int a, short int b) {
    short int a_coeffs[N_COEFFS]   = {0};
    short int b_coeffs[N_COEFFS]   = {0};
    short int res_coeffs[N_COEFFS] = {0};

    label2config(a,a_coeffs,N_COEFFS);
    label2config(b,b_coeffs,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        res_coeffs[i] = a_coeffs[i] + b_coeffs[i];
    }
    return config2label(res_coeffs,N_COEFFS);
}



void compute_sum_table(short int sum_table[][N_CONFIGS]) {
    for (int a = 0; a < N_CONFIGS; ++a) {
        for (int b = 0; b < N_CONFIGS; ++b) {
            if (a == ZERO_LABEL)      sum_table[a][b] = b;
            else if (b == ZERO_LABEL) sum_table[a][b] = a;
            else {
                sum_table[a][b] = sum_two_vectors(a,b);
            }
        }
    }
}



__attribute__((unused))
short int zero_label() {
    const short int coeffs[N_COEFFS] = {0};
    return config2label(coeffs, N_COEFFS);
}



bool is_fundamental(short int label){
    // Vector is not fundamental if k is present; this is the case if
    // label >= N_CONFIGS/2
    if (label >= N_CONFIGS/2) return false;

    short int coeffs[N_COEFFS] = {0};
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
bool unique_elements(const short int array[], short int length, short int skip) {
    for (int i = 0; i < length; ++i) {
        short int val = array[i];
        if (val == skip) continue;
        for (int j = i + 1; j < length; ++j) {
            if (array[j] == val) return false;
        }
    }
    return true;
}



__attribute__((unused))
void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            printf("%g\t",gsl_matrix_get(m,i,j));
        }
        printf("\n");
    }
}
