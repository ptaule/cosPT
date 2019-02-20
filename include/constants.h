/*
   constants.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_matrix.h>

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


#endif /* ifndef CONSTANTS_H */
