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
#include "../include/version.h"
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
        const char* x_grid_file,  /* in, file with x-grid                             */
        const char* y_grid_file,  /* in, file with y-grid                             */
        const char* data_file,    /* in, file with data (z-points)                    */
        gsl_interp_accel** x_acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_interp_accel** y_acc, /* out, gsl_interpolation accelerated lookup object */
        gsl_spline2d** spline     /* out, gsl_spline_2d of values read from file      */
        )
{
    double* x = (double*)calloc(MAX_RESOLUTION, sizeof(double));
    double* y = (double*)calloc(MAX_RESOLUTION, sizeof(double));
    double* z = (double*)calloc(MAX_RESOLUTION * MAX_RESOLUTION, sizeof(double));

    size_t x_idx  = 0;
    size_t x_size = 0;
    size_t y_idx  = 0;
    size_t y_size = 0;

    FILE* fp;
    char* line = NULL;
    size_t n = 0; // Buffer size, changed by getline()
    ssize_t read; // getline() success/error flag

    // Read x-grid
    fp = fopen(x_grid_file,"r");
    if (fp == NULL) {
        error_verbose("Could not open %s. Exiting.", x_grid_file);
    }

    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (x_idx == MAX_RESOLUTION) {
            fclose(fp);
            error_verbose("Number of points in input file %s exceeds "
                    "MAX_RESOLUTION = %d. Exiting.", x_grid_file ,MAX_RESOLUTION);
        }

        int items = sscanf(line,"%lg",&x[x_idx]);

        if (items != 1) {
            fclose(fp);
            error_verbose("Reading %s: Found row where the number of items "
                    "is not equal to one. Exiting.", x_grid_file);
        }
        x_idx++;
    }
    x_size = x_idx;
    x_idx = 0;
    fclose(fp);

    // Read y-grid
    line = NULL;
    n = 0;

    fp = fopen(y_grid_file,"r");
    if (fp == NULL) {
        error_verbose("Could not open %s. Exiting.", y_grid_file);
    }

    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (y_idx == MAX_RESOLUTION) {
            fclose(fp);
            error_verbose("Number of points in input file %s exceeds "
                    "MAX_RESOLUTION = %d. Exiting.", y_grid_file, MAX_RESOLUTION);
        }

        int items = sscanf(line,"%lg",&y[y_idx]);

        if (items != 1) {
            fclose(fp);
            error_verbose("Reading %s: Found row where the number of items "
                    "is not equal to one. Exiting.", y_grid_file);
        }
        y_idx++;
    }
    y_size = y_idx;
    y_idx = 0;
    fclose(fp);

    // Read data (z-values). Grid should match x- and y-grid
    line = NULL;
    n = 0;

    fp = fopen(data_file, "r");
    if (fp == NULL) {
        error_verbose("Could not open %s. Exiting.", data_file);
    }

    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (x_idx >= x_size) {
            free(x); free(y); free(z); free(line); fclose(fp);
            error_verbose("Dimension mismatch between x-grid and first \
                    dimension of z-grid in %s. Exiting.", data_file);
        }

        y_idx = 0;
        size_t z_idx = x_idx;
        int offset = 0;
        char* data = line;
        while (sscanf(data, " %lf%n", &z[z_idx], &offset) == 1) {
            data += offset;

            if (y_idx >= y_size) {
                free(x); free(y); free(z); free(line); fclose(fp);
            error_verbose("Dimension mismatch between y-grid and first \
                    dimension of z-grid in %s. Exiting.", data_file);
            }
            y_idx++;
            /* 1d index conversion in accordance with GSL 2d spline lib */
            z_idx = y_idx*x_size + x_idx;
        }
        x_idx++;
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



double get_wavenumber(
        const char* filename, /* in, name of file to be read                 */
        int a                 /* wavenumber (zero-based) index;              */
                              /* which index (linenumber) in list to be read */
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

    int i = 0;

    while ((read = getline(&line,&n,fp) != -1)) {
        if (!strip_line(line)) continue;

        if (i == a) {
            double wavenumber = 0;
            sscanf(line,"%lg",&wavenumber);
            fclose(fp);
            free(line);

            return wavenumber;
        }

        i++;
    }

    fclose(fp);
    free(line);

    warning_verbose("Wavenumber index given (%d) is outside range of list in "
            "%s. Using k=1.0", a, filename);
    return 1.0;
}



void write_PS(
        const char* filename,
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
    fprintf(fp,"# for k=%e (h/Mpc)\n#\n", output->k);
    fprintf(fp,"# Description: %s\n", output->description);
    fprintf(fp,"# Git hash:    %s\n", build_git_sha);
    fprintf(fp,"# Build time:  %s\n", build_git_time);

    fprintf(fp,"#\n# Settings/constants used:\n#\n");
    fprintf(fp,"# component_A, component_B              = %d, %d\n#\n",
            COMPONENT_A, COMPONENT_B);
    fprintf(fp,"# Input PS read from                    = %s\n",
            output->input_ps_file);
    fprintf(fp,"# Zeta values read from                 = %s\n",
            output->input_zeta_file);
    fprintf(fp,"# Redshift values read from             = %s\n",
            output->input_redshift_file);
    fprintf(fp,"# Exp of grow. mode eig. vals from      = %s\n",
            output->input_omega_eigvals_file);
    for (int i = 0; i < COMPONENTS; ++i) {
        fprintf(fp,"# Initial F1_%d kernel read from         = %s\n",
                i, output->ic_F1_files[i]);
    }

    fprintf(fp,"#\n");
    fprintf(fp,"# Integration limits                    = [%e,%e]\n", Q_MIN, Q_MAX);
    fprintf(fp,"# Initial/final times                   = [%e,%e]\n", ETA_I, ETA_F);
    fprintf(fp,"# Asymptotic time for initialization    = %e\n#\n", ETA_ASYMP);

    fprintf(fp,"# Neutrino mass                         = %f\n", M_NU);
    fprintf(fp,"# Neutrino fraction                     = %f\n", F_NU);
    fprintf(fp,"# Square root of Omega matter at z = 0  = %f\n#\n", SQRT_OMEGA_M);

    fprintf(fp,"# Number of time steps (total)          = %d\n", TIME_STEPS);
    fprintf(fp,"# Number of time steps (before ETA_I)   = %d\n", PRE_TIME_STEPS);
    fprintf(fp,"# Monte Carlo abstol, reltol            = %.2e, %.2e\n",
            output->cuba_epsabs, output->cuba_epsrel);
    fprintf(fp,"# Monte Carlo max num. of evals         = %.2e\n",
            output->cuba_maxevals);
    fprintf(fp,"# ODE initial step size, abstol, reltol = %.2e, %.2e, %.2e\n",
            ODE_HSTART, ODE_ATOL, ODE_RTOL);
    fprintf(fp,"# ODE routine                           = %s\n",
            TOSTRING(ODE_ROUTINE));

    /* A column consists of 12 characters; 4 whitespaces in between each column */
    fprintf(fp,"#\n#%3s%-14s","","k (h/Mpc)");

    char* corr_strings[INTEGRAND_COMPONENTS] = {"<AA>","<AB>","<BB>"};

    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) {
        fprintf(fp,"%2sP_lin %-10s", "", corr_strings[i]);
        fprintf(fp,"P_%dloop %-8s", LOOPS, corr_strings[i]);
        fprintf(fp,"err_%dloop %s", LOOPS, corr_strings[i]);
    }

    for (int i = 0; i < COMPONENTS; ++i) {
        fprintf(fp,"%2sF%d(ETA_I)%5s", "", i, "");
    }

    fprintf(fp,"\n%3s% .6e", "", output->k);

    for (int i = 0; i < INTEGRAND_COMPONENTS; ++i) {
        fprintf(fp,"%3s% .6e%3s% .6e%3s% .6e", "", output->lin_ps[i], "",
                output->non_lin_ps[i], "", output->error[i]);
    }
    for (int i = 0; i < COMPONENTS; ++i) {
        fprintf(fp,"%3s% .6e", "", output->F1_eta_i[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
}
