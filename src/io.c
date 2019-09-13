/*
   io.c

   Created by Petter Taule on 25.02.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "../include/constants.h"
#include "../include/io.h"

// Maximum input resolution
#define MAX_RESOLUTION 500


// Returns false if line is empty or starts with a '#'
static bool strip_line(char line[]) {
    char * p = line;
    size_t len = strlen(line);

    // Strip newline or carriage return
    while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
        line[--len] = 0;

    if (len == 0) return false;

    // Advance to first non-whitespace
    while (isspace(*p)) p++;

    // Skip lines beginning with #
    if (*p == '#') return false;

    return true;
}



void read_and_interpolate(
        const char* filename,   /* in, name of file to be read                      */
        gsl_interp_accel** acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_spline** spline     /* out, gsl_spline of values read from file         */
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

    double* x = (double*)calloc(MAX_RESOLUTION, sizeof(double));
    double* y = (double*)calloc(MAX_RESOLUTION, sizeof(double));

    int i = 0;
    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (i == MAX_RESOLUTION) {
            fclose(fp);
            error_verbose("Number of points in input file %s exceeds "
                    "MAX_RESOLUTION = %d. Exiting.", filename,MAX_RESOLUTION);
        }

        int items = sscanf(line,"%lg" "\t" "%lg",&x[i],&y[i]);

        if (items != 2) {
            fclose(fp);
            error_verbose("Reading %s: Found row where the number of items "
                    "is not equal to two. Exiting.", filename);
        }
        i++;
    }

    // Interpolate values
    *acc = gsl_interp_accel_alloc();
    *spline = gsl_spline_alloc(INTERPOL_TYPE, i);
    gsl_spline_init(*spline,x,y,i);

    fclose(fp);
    free(line);
    free(x);
    free(y);
}



void write_PS(
        const char* filename,
        int n_points,
        const double wavenumbers[],
        const double power_spectrum[],
        const double errors[]
        )
{
    FILE* fp;

    fp = fopen(filename,"w");
    if (fp == NULL) {
        warning_verbose("Could not open %s for writing.",filename);
    }

    fprintf(fp,"# Matter power spectrum P(k) at %d-loop (%d "
        "components)\n",LOOPS, COMPONENTS);
    fprintf(fp,"# for k=%e to %e (h/Mpc)\n",
            wavenumbers[0], wavenumbers[n_points-1]);
    fprintf(fp,"# Number of wavenumbers: %d\n", n_points);

    fprintf(fp,"#\n# Settings/constants used:\n");
    fprintf(fp,"# Number of time steps                  = %d\n",TIME_STEPS);

    fprintf(fp,"# Monte Carlo abstol, reltol            = %.2e, %.2e\n", CUBA_EPSABS, CUBA_EPSREL);
    fprintf(fp,"# Monte Carlo max num. of evals         = %.2e\n", CUBA_MAXEVAL);
    fprintf(fp,"# ODE initial step size, abstol, reltol = %.2e, %.2e, %.2e\n",
            ODE_HSTART, ODE_ATOL, ODE_RTOL);
    fprintf(fp,"# ODE routine                           = %s\n", TOSTRING(ODE_ROUTINE));
    fprintf(fp,"#\n# \tk\t\t\t\tP(k)\t\t\terror\n");

    for (int i = 0; i < n_points; ++i) {
        fprintf(fp,"\t%e\t%e\t%e\n", wavenumbers[i], power_spectrum[i], errors[i]);
    }

    fclose(fp);
}
