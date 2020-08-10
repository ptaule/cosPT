/*
   main.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include <cvode/cvode.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>  // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "include/constants.h"
#include "include/wavenumbers.h"
#include "include/tables.h"
#include "include/io.h"

typedef struct {
    short int n;
    double k;
    gsl_spline** rhs_splines;
    gsl_interp_accel** rhs_accs;
    const evolution_params_t* parameters;
} ode_input_t;


static int kernel_gradient(realtype eta, N_Vector y_vec, N_Vector dy_vec, void *ode_input);
static int jacobian(realtype eta, N_Vector y, N_Vector fy, SUNMatrix J, void
        *ode_input, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int check_flag(void *flagvalue, const char *funcname, int opt);

// This macro gives access to the individual components of the data array of an
// N Vector.
#define NV_Ith_S(v,i) ( NV_DATA_S(v)[i] )


int main (int argc, char* argv[]) {
    // Input files
    const char* zeta_file            = CLASS_PATH "zeta_of_etaD.dat";
    const char* redshift_file        = CLASS_PATH "redshift_of_etaD.dat";
    const char* wavenumber_grid_file = "input/wavenumbers_bao_zoom.dat";

#if SOUND_SPEED == CG2
    const char* effcs2_etaD_grid_file = "";
    const char* effcs2_k_grid_file    = "";
    const char* effcs2_file           = "";
    const char* omega_eigvals_file =
        "input/m_nu_" M_NU_STRING "/cg2/growing_mode_eigenvalues_etaD_ini.dat";
    const char* ic_F1_files[COMPONENTS];
    ic_F1_files[0] = "input/m_nu_" M_NU_STRING "/cg2/F1_growing_mode_etaD_-10.dat";
    ic_F1_files[1] = "input/m_nu_" M_NU_STRING "/cg2/F2_growing_mode_etaD_-10.dat";
    ic_F1_files[2] = "input/m_nu_" M_NU_STRING "/cg2/F3_growing_mode_etaD_-10.dat";
    ic_F1_files[3] = "input/m_nu_" M_NU_STRING "/cg2/F4_growing_mode_etaD_-10.dat";
#elif SOUND_SPEED == EFFCS2
    const char* effcs2_etaD_grid_file = EFFCS2_PATH "m_nu_" M_NU_STRING
        "/etaD_grid.dat";
    const char* effcs2_k_grid_file    = EFFCS2_PATH "k_grid.dat";
    const char* effcs2_file           = EFFCS2_PATH "m_nu_" M_NU_STRING
        "/effcs2_exact.dat";

    const char* omega_eigvals_file =
        "input/m_nu_" M_NU_STRING "/effcs2_exact/growing_mode_eigenvalues_etaD_ini.dat";
    const char* ic_F1_files[COMPONENTS];
    ic_F1_files[0] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F1_growing_mode_etaD_-10.dat";
    ic_F1_files[1] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F2_growing_mode_etaD_-10.dat";
    ic_F1_files[2] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F3_growing_mode_etaD_-10.dat";
    ic_F1_files[3] = "input/m_nu_" M_NU_STRING "/effcs2_exact/F4_growing_mode_etaD_-10.dat";
#endif

    // Initialize time steps in eta
    double eta[TIME_STEPS];
    initialize_timesteps(eta, ETA_I, ETA_F, ETA_ASYMP);

    evolution_params_t params = {
        .zeta_acc             = NULL,
        .zeta_spline          = NULL,
        .redshift_acc         = NULL,
        .redshift_spline      = NULL,
        .omega_eigvals_acc    = NULL,
        .omega_eigvals_spline = NULL,
        .ic_F1_accs           = {NULL},
        .ic_F1_splines        = {NULL},
        .effcs2_x_acc         = NULL,
        .effcs2_y_acc         = NULL,
        .effcs2_spline        = NULL
    };

    // Read input files and interpolate
    read_and_interpolate(redshift_file, &params.redshift_acc,
            &params.redshift_spline);
    read_and_interpolate(zeta_file, &params.zeta_acc, &params.zeta_spline);
    read_and_interpolate(omega_eigvals_file, &params.omega_eigvals_acc,
            &params.omega_eigvals_spline);

#if SOUND_SPEED == EFFCS2
    read_and_interpolate_2d(effcs2_etaD_grid_file, effcs2_k_grid_file,
            effcs2_file, &params.effcs2_x_acc, &params.effcs2_y_acc,
            &params.effcs2_spline);
#endif

    for (int i = 0; i < COMPONENTS; ++i) {
        read_and_interpolate(ic_F1_files[i], &params.ic_F1_accs[i],
                &params.ic_F1_splines[i]);
    }

    // Set up ODE input and system
    ode_input_t ode_input = {
        .n = 1,
        .k = 0.0,
        .parameters = &params,
        .rhs_splines = NULL,
        .rhs_accs = NULL,
    };

    int flag; // For checking if functions have run properly
    realtype abstol = 1e-8; // absolute tolerance of system
    realtype reltol = 1e-4; // relative tolerance of system

    // 2. Defining the length of the problem.
    sunindextype N = 4;
    // 3. Set vector of initial values.
    N_Vector y; // Problem vector.
    y = N_VNew_Serial(N);
    // 4. Create CVODE Object.
    void *cvode_mem = NULL; // Problem dedicated memory.
    SUNMatrix A;
    SUNLinearSolver LS;
    cvode_mem = CVodeCreate(CV_BDF);

    double** kernels = (double**)calloc(TIME_STEPS, sizeof(double*));
    for (int i = 0; i < TIME_STEPS; ++i) {
        kernels[i] = (double*)calloc(COMPONENTS, sizeof(double));
    }

    double k = 0;
    for (int l = 0; l < NUM_WAVENUMBERS; ++l) {
        k = wavenumbers[l];
        ode_input.k = k;

        for (int i = 0; i < COMPONENTS; ++i) {
            kernels[0][i] = gsl_spline_eval(params.ic_F1_splines[i], k,
                    params.ic_F1_accs[i]);
        }

        for (int i = 1; i < PRE_TIME_STEPS + 1; ++i) {
            for (int j = 0; j < COMPONENTS; ++j) {
                kernels[i][j] = kernels[0][j] *
                    exp(gsl_spline_eval(params.omega_eigvals_spline, k,
                                params.omega_eigvals_acc) * (eta[i] - eta[0]));
            }
        }

        NV_Ith_S(y, 0) = kernels[PRE_TIME_STEPS][0];
        NV_Ith_S(y, 1) = kernels[PRE_TIME_STEPS][1];
        NV_Ith_S(y, 2) = kernels[PRE_TIME_STEPS][2];
        NV_Ith_S(y, 3) = kernels[PRE_TIME_STEPS][3];

        // 5. Initialize CVODE solver.
        realtype eta0 = eta[PRE_TIME_STEPS]; // Initiale value of time.
        flag = CVodeInit(cvode_mem, kernel_gradient, eta0, y);
#if DEBUG==1
        if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);
#endif

        // 6. Specify integration tolerances.
        flag = CVodeSStolerances(cvode_mem, reltol, abstol);
#if DEBUG==1
        if (check_flag(&flag, "CVodeSStolerances", 1)) return(1);
#endif

        // 7. Optional inputs
        flag = CVodeSetUserData(cvode_mem, &ode_input);
#if DEBUG==1
        if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);
#endif
        flag = CVodeSetStopTime(cvode_mem, ETA_F);
#if DEBUG==1
        if (check_flag(&flag, "CVodeSetStopTime", 1)) return(1);
#endif

        /* Create dense SUNMatrix for use in linear solves */
        A = SUNDenseMatrix(N, N);
#if DEBUG==1
        if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);
#endif

        /* Create dense SUNLinearSolver object for use by CVode */
        LS = SUNLinSol_Dense(y, A);
#if DEBUG==1
        if(check_flag((void *)LS, "SUNLinSol_Dense", 0)) return(1);
#endif

        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        int retval = CVodeSetLinearSolver(cvode_mem, LS, A);
#if DEBUG==1
        if(check_flag(&retval, "CVodeSetLinearSolver", 1)) return(1);
#endif

        /* Set the user-supplied Jacobian routine Jac */
        retval = CVodeSetJacFn(cvode_mem, jacobian);
#if DEBUG==1
        if(check_flag(&retval, "CVodeSetJacFn", 1)) return(1);
#endif


        // 14. Advance solution in time.
        realtype eta_reached = 0;
        for (int i = PRE_TIME_STEPS + 1; i < TIME_STEPS; i++) {
            double eta_current = eta[i];
            flag = CVode(cvode_mem, eta_current, y, &eta_reached, CV_NORMAL);
            /* N_VPrint_Serial(y); */
#if DEBUG==1
            if(check_flag(&flag, "CVode", 1)) break;
#endif
        }
        printf("%e\t%e\t%e\t%e\t%e\n", k, NV_Ith_S(y,0), NV_Ith_S(y,1), NV_Ith_S(y,2), NV_Ith_S(y,3));
    }

    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    for (int i = 0; i < TIME_STEPS; ++i) {
        free(kernels[i]);
    }
    free(kernels);

    gsl_spline_free(params.redshift_spline);
    gsl_interp_accel_free(params.redshift_acc);
    gsl_spline_free(params.zeta_spline);
    gsl_interp_accel_free(params.zeta_acc);
    gsl_spline_free(params.omega_eigvals_spline);
    gsl_interp_accel_free(params.omega_eigvals_acc);

    for (int i = 0; i < COMPONENTS; ++i) {
        gsl_interp_accel_free(params.ic_F1_accs[i]);
        gsl_spline_free(params.ic_F1_splines[i]);
    }

    return 0;
}



static int jacobian(realtype eta, N_Vector y, N_Vector fy, SUNMatrix J, void
        *ode_input, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    ode_input_t input = *(ode_input_t*)ode_input;
    short int n = input.n;
    double k = input.k;
    const evolution_params_t* params = input.parameters;

    // Between eta_asymp and eta_i, set eta = eta_i
    if (eta < ETA_I) {
        eta = ETA_I;
    }

    // zeta always enters with prefactor 1.5, hence we redefine and multiply once
    double zeta = 1.5 * gsl_spline_eval(params->zeta_spline, eta, params->zeta_acc);

    // cs2_factor = 2/3 * 1/(omegaM a^2 H^2) * (1 + z)
#define cs2_factor 2.0/3.0 * 3e3 * 3e3 / OMEGA_M_0

    double cs2 = 0.0;
#if SOUND_SPEED == CG2
#define T_nu0 1.67734976e-4 /* Neutrino temperature today [eV] */
    // Adiabatic sound speed
    // 5/9 Zeta(5)/Zeta(3) = 7.188..
    cs2 = 7.188565369 * cs2_factor * T_nu0 * T_nu0 / (M_NU * M_NU)
        * (1 + gsl_spline_eval(params->redshift_spline, eta, params->redshift_acc));
#undef T_nu0

#elif SOUND_SPEED == EFFCS2
    // Effective sound speed
    cs2 = cs2_factor * gsl_spline2d_eval(params->effcs2_spline, eta, k,
            params->effcs2_x_acc, params->effcs2_y_acc) *
        pow( 1 + gsl_spline_eval(params->redshift_spline, eta,
                    params->redshift_acc), -1);
#endif

    SM_ELEMENT_D(J, 0, 0) = -n;
    SM_ELEMENT_D(J, 0, 1) = 1;
    SM_ELEMENT_D(J, 0, 2) = 0;
    SM_ELEMENT_D(J, 0, 3) = 0;

    SM_ELEMENT_D(J, 1, 0) = zeta * (1 - F_NU);
    SM_ELEMENT_D(J, 1, 1) = - zeta + 1 - n;
    SM_ELEMENT_D(J, 1, 2) = zeta * F_NU;
    SM_ELEMENT_D(J, 1, 3) = 0;

    SM_ELEMENT_D(J, 0, 0) = 0;
    SM_ELEMENT_D(J, 0, 1) = 0;
    SM_ELEMENT_D(J, 0, 2) = -n;
    SM_ELEMENT_D(J, 0, 3) = 1;

    SM_ELEMENT_D(J, 0, 0) = zeta * (1 - F_NU);
    SM_ELEMENT_D(J, 0, 1) = 0;
    SM_ELEMENT_D(J, 0, 2) = zeta * (F_NU - k*k*cs2);
    SM_ELEMENT_D(J, 0, 3) = - zeta + 1 - n;

#undef cs2_factor
    return 0;
}



static int kernel_gradient(realtype eta, N_Vector y_vec, N_Vector dy_vec, void *ode_input) {
    ode_input_t input = *(ode_input_t*)ode_input;
    short int n = input.n;
    double k = input.k;
    const evolution_params_t* params = input.parameters;

    realtype* y  = N_VGetArrayPointer(y_vec);
    realtype* dy = N_VGetArrayPointer(dy_vec);

    // Between eta_asymp and eta_i, set eta = eta_i
    if (eta < ETA_I) {
        eta = ETA_I;
    }

    // Interpolate RHS
    double rhs[COMPONENTS] = {0};
    // If n == 1, rhs = 0
    if (n > 1) {
        for (int i = 0; i < COMPONENTS; ++i) {
            rhs[i] = gsl_spline_eval(input.rhs_splines[i], eta, input.rhs_accs[i]);
        }
    }

    // zeta always enters with prefactor 1.5, hence we redefine and multiply once
    double zeta = 1.5 * gsl_spline_eval(params->zeta_spline, eta, params->zeta_acc);

    // cs2_factor = 2/3 * 1/(omegaM a^2 H^2) * (1 + z)
#define cs2_factor 2.0/3.0 * 3e3 * 3e3 / OMEGA_M_0

    double cs2 = 0.0;
#if SOUND_SPEED == CG2
#define T_nu0 1.67734976e-4 /* Neutrino temperature today [eV] */
    // Adiabatic sound speed
    // 5/9 Zeta(5)/Zeta(3) = 7.188..
    cs2 = 7.188565369 * cs2_factor * T_nu0 * T_nu0 / (M_NU * M_NU)
        * (1 + gsl_spline_eval(params->redshift_spline, eta, params->redshift_acc));
#undef T_nu0

#elif SOUND_SPEED == EFFCS2
    // Effective sound speed
    cs2 = cs2_factor * gsl_spline2d_eval(params->effcs2_spline, eta, k,
            params->effcs2_x_acc, params->effcs2_y_acc) *
        pow( 1 + gsl_spline_eval(params->redshift_spline, eta,
                    params->redshift_acc), -1);
#endif

    // etaD parametrization with zeta(etaD) from CLASS
    dy[0] = rhs[0] - n * y[0] + y[1];
    dy[1] = rhs[1] + zeta * (1 - F_NU) * y[0] + (- zeta + 1 - n) * y[1] + zeta * F_NU * y[2];
    dy[2] = rhs[2] - n * y[2] + y[3];
    dy[3] = rhs[3] + zeta * (1 - F_NU) * y[0] + zeta * (F_NU - k*k*cs2) * y[2] +
        (- zeta + 1 - n) * y[3];

#undef cs2_factor
    return 0;
}



// check_flag function is from the cvDiurnals_ky.c example from the CVODE
// package.
/* Check function return value...
   opt == 0 means SUNDIALS function allocates memory so check if
   returned NULL pointer
   opt == 1 means SUNDIALS function returns a flag so check if
   flag >= 0
   opt == 2 means function allocates memory so check if returned
   NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt) {
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return(1); }}

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}
