/*
   utilities.h

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stddef.h>

// Constants:
#define PI 3.14159265359
#define TWOPI 6.28318530718

#define LOOPS 1
#define N_COEFFS LOOPS+1
#define N_CONFIGS (int)(2 * pow(3,LOOPS))
#define MAX_KERNEL_NUM 2*LOOPS + 1

#define COMPONENTS 2

#define DEBUG 1

// Macros:

// debug-print is optimized away if DEBUG==false
#define debug_print(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)


// Utility functions

void label2config(short int label, short int config[], size_t size);
short int config2label(const short int config[], size_t size);

short int sum_two_vectors(short int label_a, short int label_b);
short int sum_vectors(const short int labels[], size_t size);

short int numberOfKernels(short int n, short int n_configs);


struct gsl_matrix;
void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width);

#endif /* ifndef UTILITIES_H */
