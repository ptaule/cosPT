/*
   ir_resum.cpp

   Created by Petter Taule on 15.02.2023
   Copyright (c) 2023 Petter Taule. All rights reserved.
*/

#include <algorithm>
#include <cstddef>
#include <exception>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "../include/interpolation.hpp"
#include "../include/ir_resum.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;

struct IRDampingIntParams {
    const Interpolation1D& ps_nw;
    double k_osc;

    IRDampingIntParams(
            const Interpolation1D& ps_nw,
            double k_osc
            ) : ps_nw(ps_nw), k_osc(k_osc) {}
};



double ir_Sigma_integrand(double q, void* parameters) {
    IRDampingIntParams params =
        *static_cast<IRDampingIntParams*>(parameters);
    double r = q/params.k_osc;

    return params.ps_nw(q)
        * (1 - gsl_sf_bessel_j0(r) + 2 * gsl_sf_bessel_j2(r));
}



double ir_delta_Sigma_integrand(double q, void* parameters) {
    IRDampingIntParams params =
        *static_cast<IRDampingIntParams*>(parameters);
    double r = q/params.k_osc;

    return params.ps_nw(q) * gsl_sf_bessel_j2(r);
}



void compute_ir_damping(
        const Interpolation1D& ps_nw,
        double k_min,
        double k_s,
        double k_osc,
        double& Sigma2,        /* out */
        double& delta_Sigma2,  /* out */
        size_t sub_regions,
        double atol,
        double rtol,
        int key
        )
{
    gsl_integration_workspace* workspace;
    workspace = gsl_integration_workspace_alloc(sub_regions);

    IRDampingIntParams params(ps_nw, k_osc);

    gsl_function F;
    F.function = ir_Sigma_integrand;
    F.params = static_cast<void*>(&params);

    double error;
    int status = gsl_integration_qag(&F, k_min, k_s, atol, rtol, sub_regions,
            key, workspace, &Sigma2, &error);
    Sigma2 *= 4.0 * PI / 3.0;

    if (status != 0) {
        throw std::runtime_error("IR damping function Sigma^2 integration \
                failed with error code" + std::to_string(status));
    }

    F.function = ir_delta_Sigma_integrand;
    status = gsl_integration_qag(&F, k_min, k_s, atol, rtol, sub_regions,
            key, workspace, &delta_Sigma2, &error);
    delta_Sigma2 *= 4.0 * PI;

    if (status != 0) {
        throw std::runtime_error("IR damping function deltaSigma^2 \
                integration failed with error code"
                + std::to_string(status));
    }

    gsl_integration_workspace_free(workspace);
}
