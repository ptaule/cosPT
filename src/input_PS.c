/*
   input_PS.c

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "../include/constants.h"
#include "../include/input_PS.h"



void read_input_PS(
        const char* filename,   /* in, power spectrum file                          */
        gsl_interp_accel** acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_spline** spline     /* out, copy of data read from file                 */
        )
{
    FILE* fp;
    char* line = NULL;
    size_t n = 0; // Buffer size, changed by getline()
    ssize_t read; // getline() success/error flag

    fp = fopen(filename,"r");
    if (fp == NULL) {
        error_verbose("Could not open %s. Exiting.",filename);
    }

    double* wavenumber     = (double*)malloc(sizeof(double) * MAX_RESOLUTION);
    double* power_spectrum = (double*)malloc(sizeof(double) * MAX_RESOLUTION);

    int i = 0;

    // Read line-by-line
    while ((read = getline(&line,&n,fp) != -1)) {
        char * p = line;
        size_t len = strlen(line);

        // Strip newline or carriage return
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = 0;

        if (len == 0) continue;

        // Advance to first whitespace
        while (isspace(*p)) p++;

        // Skip lines beginning with #
        if (*p == '#') continue;

        if(i == MAX_RESOLUTION) {
            error_verbose("Number of points in input file %s exceeds "
                    "MAX_RESOLUTION = %d. Exiting.", filename,MAX_RESOLUTION);
        }

        int items = sscanf(line,"%lg" "\t" "%lg",&wavenumber[i],&power_spectrum[i]);

        if (items != 2) {
            error_verbose("Reading %s: Found row where the number of items "
                    "does not equal two. Exiting.", filename);
        }

        i++;
    }

    debug_print("RESOLUTION    = %d\n",i);

    // Interpolate values
    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(INTERPOL_TYPE, i);
    gsl_spline_init(*spline,wavenumber,power_spectrum,i);

    fclose(fp);
    free(line);
    free(wavenumber);
    free(power_spectrum);
}
