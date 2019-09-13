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

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

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



void read_and_interpolate_2d(
        const char* filename,     /* in, name of file to be read                      */
        gsl_interp_accel** x_acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_interp_accel** y_acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_spline2d** spline     /* out, gsl_spline_2d of values read from file      */
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
    double* z = (double*)calloc(MAX_RESOLUTION * MAX_RESOLUTION, sizeof(double));

    size_t x_size = 0;
    size_t y_size = 0;
    size_t i = 0;

    // Read first two (non-empty, non-comment) lines
    while ((read = getline(&line,&n,fp) != -1) && i < 2) {
        if (!strip_line(line)) continue;

        /* First non-empty, non-comment line is the x grid */
        if (i == 0) {
            size_t j = 0;
            int offset = 0;
            char* data = line;
            while (sscanf(data, " %lf%n", &x[j], &offset) == 1) {
                j++;
                data += offset;

                if (j >= MAX_RESOLUTION) {
                    free(x); free(y); free(z); free(line); fclose(fp);
                    error_verbose("Number of grid points defined in first "
                            "(non-empty, non-comment) line in %s exceeds "
                            "MAX_RESOLUTION = %d. " "Exiting.",
                            filename,MAX_RESOLUTION);
                }
            }
            x_size = j;

        }
        /* Second non-empty, non-comment line is the y grid */
        if (i == 1) {
            size_t j = 0;
            int offset = 0;
            char* data = line;
            while (sscanf(data, " %lf%n", &y[j], &offset) == 1) {
                j++;
                data += offset;

                if (j >= MAX_RESOLUTION) {
                    free(x); free(y); free(z); free(line); fclose(fp);
                    error_verbose("Number of grid points defined in second "
                            "(non-empty, non-comment) line in %s exceeds "
                            "MAX_RESOLUTION = %d. " "Exiting.",
                            filename,MAX_RESOLUTION);
                }
            }
            y_size = j;
        }
        i++;
    }

    /* The remaining lines contain the z values*/
    i = 0;
    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (i >= x_size) {
            free(x); free(y); free(z); free(line); fclose(fp);
            error_verbose("x grid size and z values provided does match in "
                    "%s. Exiting.", filename);
        }
        size_t y_index = 0;
        size_t z_index = i;
        int offset = 0;
        char* data = line;
        while (sscanf(data, " %lf%n", &z[z_index], &offset) == 1) {
            data += offset;

            if (y_index >= y_size) {
                free(x); free(y); free(z); free(line); fclose(fp);
                error_verbose("y grid size and z values provided does match in "
                        " %s. Exiting.", filename);
            }
            y_index++;
            /* 1d index conversion in accordance with GSL 2d spline lib */
            z_index = y_index*x_size + i;
        }
        i++;
    }

    // Interpolate values
    *x_acc = gsl_interp_accel_alloc();
    *y_acc = gsl_interp_accel_alloc();
    *spline = gsl_spline2d_alloc(INTERPOL_2D_TYPE, x_size, y_size);
    gsl_spline2d_init(*spline, x, y, z, x_size, y_size);

    free(x);
    free(y);
    free(z);
    free(line);
    fclose(fp);
}



void write_PS(
        const char* filename,
        int n_points,
        const output_t* output
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
            output->wavenumbers[0], output->wavenumbers[n_points-1]);
    fprintf(fp,"# Number of wavenumbers: %d\n", n_points);

    fprintf(fp,"#\n# Settings/constants used:\n");
    fprintf(fp,"# Number of time steps                  = %d\n",TIME_STEPS);

    fprintf(fp,"# Monte Carlo abstol, reltol            = %.2e, %.2e\n", CUBA_EPSABS, CUBA_EPSREL);
    fprintf(fp,"# Monte Carlo max num. of evals         = %.2e\n", CUBA_MAXEVAL);
    fprintf(fp,"# ODE initial step size, abstol, reltol = %.2e, %.2e, %.2e\n",
            ODE_HSTART, ODE_ATOL, ODE_RTOL);
    fprintf(fp,"# ODE routine                           = %s\n", TOSTRING(ODE_ROUTINE));

    /* A column consists of 12 characters, 4 whitespaces in between each column */
    fprintf(fp,"#\n# %3s%20s%7s","k", "P_lin","");

    for (int i = 1; i <= LOOPS; ++i) {
        fprintf(fp,"%4sP_%dloop%5s","",i,"");
    }
    for (int i = 1; i <= LOOPS; ++i) {
        fprintf(fp,"%4serror_%dloop","",i);
        if (i != LOOPS) {
            fprintf(fp," ");
        }
    }
    fprintf(fp,"\n");

    for (int i = 0; i < n_points; ++i) {
        fprintf(fp,"%4s%e%4s%e%4s%e%4s%e\n", "",
                output->wavenumbers[i], "", output->lin_ps[i], "",
                output->non_lin_ps[i], "", output->errors[i]);
    }

    fclose(fp);
}
