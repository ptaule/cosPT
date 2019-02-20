/*
   utilities.h

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include "constants.h"

// Utility functions

void label2config(short int label, short int config[], size_t size);
short int config2label(const short int config[], size_t size);

short int zero_label();

bool is_fundamental(short int label);
bool unique_elements(const short int array[],size_t length, short int skip);

short int sum_two_vectors(short int label_a, short int label_b);
short int sum_vectors(const short int labels[], size_t size);

void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width);

#endif /* ifndef UTILITIES_H */
