/*
   power_spectrum_io.h

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef POWER_SPECTRUM_IO_H
#define POWER_SPECTRUM_IO_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

struct gsl_spline;
struct gsl_interp_accel;

void read_PS(const char* filename, gsl_interp_accel** acc, gsl_spline** spline);
void write_PS(const char* filename, int n_points, const double wavenumbers[], const double power_spectrum[]);
#endif /* ifndef POWER_SPECTRUM_IO_H */
