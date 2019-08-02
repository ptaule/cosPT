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
#define LOOPS 1
#endif

#define COMPONENTS 2

// Number of evaluation points,
#define N_POINTS 50
#define K_MIN 1e-2
#define K_MAX 1e1

// Integration limits
#define Q_MIN 1e-4
#define Q_MAX 1e2

// CUBA settings
#define CUBA_NVEC 1
#define CUBA_EPSREL 1e-3
#define CUBA_EPSABS 1e-12
#define CUBA_VERBOSE 0
#define CUBA_LAST 4
#define CUBA_SEED 0
#define CUBA_MINEVAL 0
#define CUBA_MAXEVAL 1e6

#define CUBA_MAXCORES 4

#define CUBA_STATEFILE NULL
#define CUBA_SPIN NULL

#define CUBA_NNEW 1000
#define CUBA_NMIN 2
#define CUBA_FLATNESS 25.

// Which GSL interpolation routine to use
#define INTERPOL_TYPE gsl_interp_cspline

// Maximum input power spectrum resolution
#define MAX_RESOLUTION 200

// Debug modes:
// 1: Perform additional checks during runtime.
// 2: In addition, print info during runtime.
#ifndef DEBUG
#define DEBUG 0
#endif

// Variable precision
typedef long double vfloat;
#define vfloat_fmt   "%Le"
#define matrix_t     gsl_matrix_long_double
#define matrix_alloc gsl_matrix_long_double_alloc
#define matrix_set   gsl_matrix_long_double_set
#define matrix_get   gsl_matrix_long_double_get
#define matrix_free  gsl_matrix_long_double_free

// Various colors for debug output

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// debug_print() is optimized away if DEBUG==0
#define debug_print(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
#define warning(fmt) \
    fprintf(stderr, ANSI_COLOR_BLUE "%s:%d:\tWarning: " fmt ANSI_COLOR_RESET "\n", \
            __FILE__, __LINE__);
#define warning_verbose(fmt, ...) \
    fprintf(stderr, ANSI_COLOR_BLUE "%s:%d:\tWarning: " fmt ANSI_COLOR_RESET "\n", \
            __FILE__, __LINE__, __VA_ARGS__);
#define error(fmt) \
    {fprintf(stderr, ANSI_COLOR_RED "%s:%d:\tError: " fmt ANSI_COLOR_RESET "\n", \
            __FILE__, __LINE__); exit(EXIT_FAILURE); }
#define error_verbose(fmt, ...) \
    {fprintf(stderr, ANSI_COLOR_RED "%s:%d:\tError: " fmt ANSI_COLOR_RESET "\n", \
            __FILE__, __LINE__, __VA_ARGS__); exit(EXIT_FAILURE); }

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

// Constants:
#define PI 3.14159265359
#define TWOPI 6.28318530718

// Constants defined depending on number of loops
#define N_COEFFS LOOPS+1
#define N_DIMS   3*LOOPS-1

#if LOOPS==1
#define N_CONFIGS      6
#define N_KERNEL_ARGS  3
#define N_KERNELS     16
#define ZERO_LABEL     1
#define N_DIAGRAMS     2
#endif

#if LOOPS==2
#define N_CONFIGS      18
#define N_KERNEL_ARGS   5
#define N_KERNELS     160
#define ZERO_LABEL      4
#define N_DIAGRAMS      4
#endif

#if LOOPS==3
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
