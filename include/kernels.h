/*
   kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef KERNELS_H
#define KERNELS_H

#include "constants.h"

typedef struct {
    vfloat value;
    bool computed;
} kernel_value;

struct gsl_matrix;

void compute_scalar_products(const double k, const double Q, const double mu, matrix_vfloat* scalar_products);

void compute_alpha_beta_tables(const double k, const double Q, const double mu, gsl_matrix* alpha, gsl_matrix* beta);

short int kernel_index_from_fundamental(short int argument);
void kernel_index_from_arguments(const short int arguments[], short int* index, short int* n);
short int combined_kernel_index(short int argument_index,short int component);

#endif /* ifndef KERNELS_H */
