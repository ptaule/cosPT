/*
   ir_resum.cpp

   Created by Petter Taule on 15.02.2023
   Copyright (c) 2023 Petter Taule. All rights reserved.
*/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include "../include/interpolation.hpp"
#include "../include/ir_resum.hpp"
#include "../include/utilities.hpp"

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

using std::size_t;
using std::pow;
using std::sqrt;
using std::sin;
using std::cos;


/* Discrete sine transform II */
void DST_II(Vec1D<double>& data) {
    size_t N = data.size();
    /* Check that N is larger than 0 and a power of 2 */
    if (N == 0 || ((N & (N-1)) != 0)) {
        throw std::invalid_argument("DST(): data.size() is either 0 or not a"
                                  "power of 2.");
    }

    Vec1D<double> input(data);
    /* Make twice as long by adding zeros. (See
     * https://en.wikipedia.org/wiki/Discrete_sine_transform ; to convert the
     * expression to the equivalent FFT sum, we need twice as long an input) */
    input.insert(input.end(), N, 0);

    /* FFT of purely real input */
    gsl_fft_real_radix2_transform(input.data(), 1, 2*N);
    /* Normalization factor 1/sqrt(2*N) */
    double norm = SQRT2 * sqrt(N);
    for (auto& el : input) el /= norm;

    for (size_t i = 0; i < N; ++i) {
        /* gsl_fft_real_radix2_transform uses the fact that the input is real
         * to store real/imag output at specific locations, see GSL
         * documentation */
        double re = input.at(i + 1);
        double im = -input.at(2*N - i - 1);
        double arg = static_cast<double>(i + 1) * PI /
            (2.0 * static_cast<double>(N));

        data.at(i) = SQRT2 * (cos(arg) * im + sin(arg) * re);
    }
}



/* Discrete sine transform III (inverse of DST_II) */
void DST_III(Vec1D<double>& data) {
    size_t N = data.size();
    /* Check that N is larger than 0 and a power of 2 */
    if (N == 0 || ((N & (N-1)) != 0)) {
        throw std::invalid_argument("DST(): data.size() is either 0 or not a"
                                  "power of 2.");
    }

    Vec1D<double> input(data);
    /* Change last element to 0 and make four times as long by adding zeros.
     * (See https://en.wikipedia.org/wiki/Discrete_sine_transform ; to convert
     * the expression to the equivalent FFT sum, we need twice as long an
     * input) */
    input.back() = 0;
    input.insert(input.end(), 3*N, 0);

    /* FFT of purely real input */
    gsl_fft_real_radix2_transform(input.data(), 1, 4*N);
    /* Normalization factor 1/sqrt(2*N) */
    double norm = 2 * sqrt(N);
    for (auto& el : input) el /= norm;

    /* Store (normalized) last element of data which is explicitly added in the
     * final sum. (See https://en.wikipedia.org/wiki/Discrete_sine_transform )
     * */
    double last = data.back() / sqrt(N);

    for (size_t i = 0; i < N; ++i) {
        /* gsl_fft_real_radix2_transform uses the fact that the input is real
         * to store real/imag output at specific locations, see GSL
         * documentation.
         * We need fourier index 2i+1 */
        size_t idx = 2*i + 1;
        double re = input.at(idx);
        double im = -input.at(4*N - idx);
        double arg = static_cast<double>(idx) * PI /
            (2.0 * static_cast<double>(N));

        data.at(i) = 4 * (cos(arg) * im + std::sin(arg) * re) +
            pow(-1,i) * last;
    }
}



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
