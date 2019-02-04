/*
   utilities.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf.h>

#include "utilities.h"


void label2config(
        short int label,    // in, label element in [0, N_CONFIG - 1]
        short int config[], // out, array of base 3 digits
        size_t size   // in, size of coefficients
        )
{
    // Q-coefficients are element in (-1,0,1)
    for (size_t i = 0; i < size - 1; ++i) {
        config[i] = (label % 3) - 1;
        label /= 3;
    }

    // k-coefficient (last coefficient) is element in (0,1)

    config[size - 1] = label % 3;

}

short int config2label(
        const short int config[], // in, array of base 3 digits
        size_t size         // in, size of coefficients
        )
{
    short int label = 0;

    // For Q-coefficients, add 1 to obtain range (0,1,2)
    for (size_t i = 0; i < size - 1; ++i) {
        label += (config[i] + 1) * pow(3,i);
    }

    // For k-coefficient, do _not_ add 1, range is (0,1)
    label += config[size-1] * pow(3,size - 1);

    return label;
}


short int sum_two_vectors(short int label_a, short int label_b) {

    short int a_coeffs[N_COEFFS]   = {};
    short int b_coeffs[N_COEFFS]   = {};
    short int res_coeffs[N_COEFFS] = {};

    label2config(label_a,a_coeffs,N_COEFFS);
    label2config(label_b,b_coeffs,N_COEFFS);

    for (int i = 0; i < N_COEFFS; ++i) {
        res_coeffs[i] = a_coeffs[i] + b_coeffs[i];
    }

    // If DEBUG==true, check that sum is an appropriate vector configuration, i.e. that Q-coefficients are elements of (-1,0,1) and k-coefficient is an element of (0,1)
#if DEBUG
    for (int i = 0; i < N_COEFFS - 1; ++i) {
        short int c = res_coeffs[i];
        if (!(c == -1 || c == 0 || c == 1))
            fprintf(stderr, "Warning: Sum of vectors with labels (%d,%d) does not \
correspond to an appropriate configuration.\n", label_a,label_b);
    }
    short int c = res_coeffs[N_COEFFS - 1];
    if (!(c == 0 || c == 1))
        fprintf(stderr, "Warning: Sum of vectors with labels (%d,%d) does not \
correspond to an appropriate configuration.\n", label_a,label_b);
#endif

    return config2label(res_coeffs,N_COEFFS);
}

short int sum_vectors(const short int labels[], size_t size) {
    if (size == 1) return labels[0];

    short int result = sum_two_vectors(labels[0],labels[1]);

    for (size_t i = 2; i < size; ++i) {
        result = sum_two_vectors(result,labels[i]);
    }
    return result;
}

short int numberOfKernels(short int n, short int n_configs) {
    short int number = 0;

    for (int i = 2; i <= n; ++i) {
        // Binomial coefficient (n_configs + i - 1,i) determines d.o.f. of a
        // symmetric i-tensor of size n_configs
        number += gsl_sf_choose(n_configs + i - 1, i);
    }
    return number;
}

void print_gsl_matrix(const gsl_matrix* m, size_t height, size_t width) {
    for (size_t i = 0; i < height; ++i) {
        for (size_t j = 0; j < width; ++j) {
            printf("%g\t",gsl_matrix_get(m,i,j));
        }
        printf("\n");
    }
}
