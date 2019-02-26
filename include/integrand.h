/*
   integrand.h

   Created by Petter Taule on 24.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_H
#define INTEGRAND_H


typedef struct {
    short int l;
    short int r;
    short int m;
} diagram_t;

int diagram_factor(const diagram_t* diagram);
int integrand_symmetrization_factor(const diagram_t* diagram);
vfloat compute_k1(short int m, const vfloat bare_scalar_products[][N_COEFFS]);
int heaviside_theta(short int m, vfloat k1, const vfloat Q_magnitudes[]);

#endif /* ifndef INTEGRAND_H */
