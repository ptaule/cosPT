/*
   io.h

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef IO_H
#define IO_H

#include <gsl/gsl_spline.h>

typedef struct {
    const char* input_ps_file;
    int n_points;
    double* wavenumbers;
    double* non_lin_ps;
    double* errors;
} output_t;

void read_and_interpolate(const char* filename, gsl_interp_accel** acc,
        gsl_spline** spline);
void write_PS(const char* filename, const output_t* output);

#endif /* ifndef IO_H */
