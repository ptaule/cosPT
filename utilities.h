/*
   utilities.h

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <stddef.h>
#include <stdbool.h>

// Constants:
#define PI 3.14159265359
#define TWOPI 6.28318530718

#define LOOPS 1
#define N_COEFFS (LOOPS+1)
#define N_CONFIGS (int)(2 * pow(3,LOOPS))
#define N_KERNEL_ARGS (2 * LOOPS + 1)
#define N_KERNELS (int)((pow(3,LOOPS) + 1) * pow(4,LOOPS))

#define COMPONENTS 2

#define ZERO_LABEL zero_label()

#define DEBUG true

// Variable precision
#define matrix_vfloat gsl_matrix
#define vector_vfloat gsl_vector
typedef double vfloat;

// Macros:

// debug-print is optimized away if DEBUG==false
#define debug_print(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
#define warning(fmt) \
                fprintf(stderr, "%s:%d,\tWarning: " fmt "\n", __FILE__, __LINE__);
#define warning_verbose(fmt, ...) \
                fprintf(stderr, "%s:%d,\tWarning: " fmt "\n",__FILE__, __LINE__, __VA_ARGS__);

// Utility functions

void label2config(short int label, short int config[], size_t size);
short int config2label(const short int config[], size_t size);

short int zero_label();

bool is_fundamental(short int label);
bool unique_elements(const short int array[],size_t length, short int skip);

short int sum_two_vectors(short int label_a, short int label_b);
short int sum_vectors(const short int labels[], size_t size);


struct gsl_matrix;
void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width);


typedef struct {
    vfloat value;
    bool computed;
} kernel_value;

#endif /* ifndef UTILITIES_H */
