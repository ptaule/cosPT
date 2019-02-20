/*
   spt_kernels.h

   Created by Petter Taule on 18.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef SPT_KERNELS_H
#define SPT_KERNELS_H

#include <gsl/gsl_matrix.h>

#include "constants.h"
#include "kernels.h"

vfloat partial_SPT_sum(
        const short int arguments[],
        short int component,
        const matrix_vfloat* alpha,
        const matrix_vfloat* beta,
        kernel_value* kernels,
        const short int n,
        const short int m,
        const short int a,
        const short int b
        );

vfloat compute_SPT_kernel(
        const short int arguments[],
        short int component,
        const matrix_vfloat* alpha,
        const matrix_vfloat* beta,
        kernel_value* kernels
        );


#endif /* ifndef SPT_KERNELS_H */
