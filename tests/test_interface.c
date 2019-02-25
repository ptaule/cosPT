/*
   test_interface.c

   Created by Petter Taule on 20.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdbool.h>
#include <gsl/gsl_matrix.h>

#include "../include/constants.h"

#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/spt_kernels.h"


vfloat SPTkernel1Loop(int n, int component, vfloat k, vfloat q, vfloat mu) {
    short int args[N_KERNEL_ARGS] = {};

    if (n == 1) return 1.0;

    if (n == 2) {
        args[0] = 3;   // 3 <-> k - q
        args[1] = 2;   // 2 <-> q
        args[2] = ZERO_LABEL;
    }

    if (n == 3) {
        args[0] = 4;   // 4 <-> k
        args[1] = 2;   // 2 <-> q
        args[2] = 0;   // 0 <-> -q
    }

    vfloat Q_magnitudes[N_COEFFS - 1] = {};
    vfloat cos_theta[N_COEFFS - 1]  = {};
    vfloat phi[N_COEFFS - 2];

    Q_magnitudes[0] = q;
    cos_theta[0]    = mu;

    matrix_vfloat* alpha = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);
    matrix_vfloat* beta = gsl_matrix_alloc(N_CONFIGS,N_CONFIGS);

    compute_alpha_beta_tables(k,Q_magnitudes,cos_theta,phi,alpha,beta);

    // Allocate space for kernels (calloc also initializes values to 0)
    kernel_value* kernels = (kernel_value*)calloc(COMPONENTS * N_KERNELS, sizeof(kernel_value));

    vfloat value = compute_SPT_kernel(args,component,alpha,beta,kernels);

    // Free allocated memory
    free(kernels);
    gsl_matrix_free(alpha);
    gsl_matrix_free(beta);

    return value;
}
