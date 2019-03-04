/*
   constants.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_matrix.h>

// Parameters (if not set by compile options)
#ifndef LOOPS
#define LOOPS 2
#endif
#ifndef COMPONENTS
#define COMPONENTS 2
#endif

// Which GSL interpolation routine to use
#define INTERPOL_TYPE gsl_interp_cspline

// Maximum input power spectrum resolution
#define MAX_RESOLUTION 200

// Debug mode. Performs additional checks during runtime
#ifndef DEBUG
#define DEBUG 1
#endif

// Variable precision
typedef double vfloat;
#define matrix_vfloat gsl_matrix
#define vector_vfloat gsl_vector
#define vfloat_fmt "%lg"

// Macros:

#define interpolate(k) gsl_spline_eval(spline,k,acc)

// debug-print is optimized away if DEBUG==0
#define debug_print(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
#define warning(fmt) \
                fprintf(stderr, "%s:%d,\tWarning: " fmt "\n", __FILE__, __LINE__);
#define warning_verbose(fmt, ...) \
                fprintf(stderr, "%s:%d,\tWarning: " fmt "\n",__FILE__, __LINE__, __VA_ARGS__);
#define error(fmt) \
    {fprintf(stderr, "%s:%d,\tWarning: " fmt "\n", __FILE__, __LINE__); \
        exit(EXIT_FAILURE); }
#define error_verbose(fmt, ...) \
    {fprintf(stderr, "%s:%d,\tWarning: " fmt "\n",__FILE__, __LINE__, __VA_ARGS__); \
        exit(EXIT_FAILURE); }


// Constants:
#define PI 3.14159265359
#define TWOPI 6.28318530718

// Constants defined depending on number of loops
#if LOOPS==1
#define N_COEFFS       2
#define N_CONFIGS      6
#define N_KERNEL_ARGS  3
#define N_KERNELS     16
#define ZERO_LABEL     1
#define N_DIAGRAMS     2
#endif

#if LOOPS==2
#define N_COEFFS        3
#define N_CONFIGS      18
#define N_KERNEL_ARGS   5
#define N_KERNELS     160
#define ZERO_LABEL      4
#define N_DIAGRAMS      4
#endif

#if LOOPS==3
#define N_COEFFS         4
#define N_CONFIGS       54
#define N_KERNEL_ARGS    7
#define N_KERNELS     1792
#define ZERO_LABEL      13
#define N_DIAGRAMS       6
#endif

/* General expressions:
#define N_COEFFS (LOOPS+1)
#define N_CONFIGS (int)(2 * pow(3,LOOPS))
#define N_KERNEL_ARGS (2 * LOOPS + 1)
#define N_KERNELS (int)((pow(3,LOOPS) + 1) * pow(4,LOOPS))
#define ZERO_LABEL zero_label()
*/

#endif /* ifndef CONSTANTS_H */
