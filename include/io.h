/*
   io.h

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef IO_H
#define IO_H

#include <gsl/gsl_spline.h>

void read_PS(const char* filename, gsl_interp_accel** acc, gsl_spline** spline);
void write_PS(const char* filename,
        int n_points,
        const double wavenumbers[],
        const double power_spectrum[],
        const double errors[]
        );

#endif /* ifndef IO_H */
