/*
   integrand.c

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>
#include <gsl/gsl_sf.h>

#include "../include/constants.h"

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


