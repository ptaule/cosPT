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



void print_label(short int label) {
    if (label == ZERO_LABEL) return;

    short int config[N_COEFFS];
    label2config(label,config,N_COEFFS);

    if (config[N_COEFFS - 1] == 1) printf("k");
    for (int i = 0; i < LOOPS; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) printf("-Q%d",i+1);
        else if (config[i] == 1)  printf("+Q%d",i+1);
    }
}



void print_labels(const short int labels[])
{
    printf("(");
    for (int i = 0; i < N_KERNEL_ARGS; ++i) {
        if (labels[i] == ZERO_LABEL) continue;
        print_label(labels[i]);
        printf(", ");
    }
    printf(")");
}



void print_gsl_matrix(const matrix_t* m, size_t height, size_t width) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            printf(vfloat_fmt "\t",matrix_get(m,i,j));
        }
        printf("\n");
    }
}
