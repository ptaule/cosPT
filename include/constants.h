/*
   constants.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

// Parameters (if not set by compile options)
#ifndef LOOPS
#define LOOPS 1
#endif

#define COMPONENTS 2
#define TIME_STEPS 100

// Min/max wavenumber
#define K_MIN 1e-4
#define K_MAX 1e2

// Number of evaluation points between k=K_MIN and k=K_MAX
#define N_POINTS 139

// Debug modes (if not set by compile options):
//
// 1: Perform additional checks during runtime.
// 2: In addition, print info during runtime.
#ifndef DEBUG
#define DEBUG 0
#endif

// Which GSL interpolation routine to use
#define INTERPOL_TYPE gsl_interp_cspline
// Which GSL ODE routine to use
#define ODE_ROUTINE gsl_odeiv2_step_rkf45
// GSL ODE routine parameters
#define ODE_HSTART 1e-4 /* Initial step size */
#define ODE_ATOL 1e-6
#define ODE_RTOL 1e-4

// Variable precision
typedef long double vfloat;
#define vfloat_fmt   "%Le"
#define matrix_t     gsl_matrix_long_double
#define matrix_alloc gsl_matrix_long_double_alloc
#define matrix_set   gsl_matrix_long_double_set
#define matrix_get   gsl_matrix_long_double_get
#define matrix_free  gsl_matrix_long_double_free

// Parameters type
typedef struct {
    vfloat eta_i;
    vfloat eta_f;
    gsl_interp_accel* zeta_acc;
    gsl_spline* zeta_spline;
    gsl_matrix* omega;
} evolution_params_t;

// Constants:
#define PI    3.14159265359
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

// Various colors for debug output

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

// Warning/error commands which print filename, line number and (optionally) a
// message.
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

#endif /* ifndef CONSTANTS_H */
