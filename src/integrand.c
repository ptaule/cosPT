/*
   integrand.c

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "../include/constants.h"
#include "../include/utilities.h"
#include "../include/kernels.h"
#include "../include/spt_kernels.h"



int diagram_factor(short int l, short int r, short int m) {
    int numerator = gsl_sf_fact(2*l + m) * gsl_sf_fact(2*r + m);
    int denominator = pow(2,l+r) * gsl_sf_fact(l) * gsl_sf_fact(r) * gsl_sf_fact(m);

    return numerator/denominator;
}



int integrand_symmetrization_factor(short int l, short int r, short int m) {
    int numerator = gsl_sf_fact(LOOPS) * pow(2,m-1);
    int denominator = gsl_sf_fact(m-1) * gsl_sf_fact(l) * gsl_sf_fact(r);

    return numerator/denominator;
}



vfloat compute_k1(short int m, const vfloat bare_scalar_products[][N_COEFFS]) {
    vfloat k1 = bare_scalar_products[N_COEFFS - 1][N_COEFFS - 1];
    for (int i = 2; i <= m; ++i) {
        k1 += bare_scalar_products[i-2][i-2];
        k1 -= 2 * bare_scalar_products[N_COEFFS - 1][i-2];
    }

    for (int i = 2; i <= m; ++i) {
        for (int j = 2; j < i; ++j) {
            k1 += 2 * bare_scalar_products[i-2][j-2];
        }
    }

    return k1;
}



int heaviside_theta(short int m, vfloat k1, const vfloat Q_magnitudes[]) {
    // Heaviside-theta (k1 - k2)
    if (m == 2) {
        if (k1 <= Q_magnitudes[0]) return 0;
        return 2;
    }

    // Heaviside-theta (k2 - k3) etc.
    // Note that (assuming m >= 2), k2 = Q_magnitudes[0] etc.
    for (int i = 3; i <= m; ++i) {
        if (Q_magnitudes[i-2] <= Q_magnitudes[i-1]) return 0;
    }
    return gsl_sf_fact(m);
}
