/*
   input_PS.h

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#ifndef INPUT_PS_H
#define INPUT_PS_H

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

void read_input_PS(const char* filename, gsl_interp_accel** acc, gsl_spline** spline);

#endif /* ifndef INPUT_PS_H */
