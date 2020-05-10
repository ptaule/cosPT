/*
   constants.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_spline.h>

// Path to input files
#define INPUT_PATH  "/space/ge52sir/CLASS/massive_nu_0.07eV/"

// Number of loops (if not set by compile options)
#ifndef LOOPS
#define LOOPS 1
#endif

#define COMPONENTS 2
#define SPT_COMPONENTS 2

#ifndef TIME_STEPS
#define TIME_STEPS 100
#endif

// Integration limits
#define Q_MIN 1e-4
#define Q_MAX 6.5e1

// Initial/final times (eta is typically defined as the log of a or the growth
// rate D+)
#define ETA_F 0.0

// Which GSL interpolation routine to use
#define INTERPOL_TYPE gsl_interp_cspline
#define INTERPOL_2D_TYPE gsl_interp2d_bilinear
// Which GSL ODE routine to use
#define ODE_ROUTINE gsl_odeiv2_step_rkf45
// GSL ODE routine parameters
#define ODE_HSTART 1e-3 /* Initial step size */
#define ODE_ATOL 1e-6
#define ODE_RTOL 1e-4

// Debug modes (if not set by compile options):
//
// 1: Perform additional checks during runtime.
// 2: In addition, print info during runtime.
#ifndef DEBUG
#define DEBUG 0
#endif

// Variable precision
typedef long double vfloat;
#define vfloat_fmt   "%Le"

// Parameters type
typedef struct {
    gsl_interp_accel* zeta_acc;
    gsl_spline* zeta_spline;
    gsl_interp_accel* ic_perturb_accs[1];
    gsl_spline* ic_perturb_splines[1];
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
