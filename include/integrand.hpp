/*
   integrand.hpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_HPP
#define INTEGRAND_HPP

#include <cmath>
#include <functional>
#include <string>
#include <utility>

#include "utilities.hpp"
#include "interpolation.hpp"

class IntegrandTables;
class PowerSpectrumDiagram;
class BiSpectrumDiagram;
class IRresumSettings;

class InputPowerSpectrum {
    private:
        Interpolation1D ps;
        Interpolation1D ps_nw;

        double q_min = 0;
        double q_max = 0;

        double ir_resum = false;
        double rsd = false;
        double rsd_growth_f = 0;
        double Sigma2 = 0;
        double delta_Sigma2 = 0;

        std::function<double(double)> Sigma2_tot;
        std::function<double(double, double)> evaluate;

        double Sigma2_total_rsd(double mu) const {
            double mu2 = mu*mu;
            return (1 + rsd_growth_f * mu2 * (2 + rsd_growth_f)) * Sigma2 +
                rsd_growth_f * rsd_growth_f * mu2 * (mu2 - 1) * delta_Sigma2;
        };

        /* _n0 to _n2 functions differ due to different correction of overcounting
         * of IR contributions at different PT and loop orders */
        double IR_resum_n0(double q, double mu) const {
            double ps_nw_q = ps_nw(q, q_min, q_max);
            double ps_w_q = ps(q, q_min, q_max) - ps_nw_q;
            double k2S2 = q * q * Sigma2_tot(mu);
            return ps_nw_q + std::exp(-k2S2) * ps_w_q;
        }
        double IR_resum_n1(double q, double mu) const {
            double ps_nw_q = ps_nw(q, q_min, q_max);
            double ps_w_q = ps(q, q_min, q_max) - ps_nw_q;
            double k2S2 = q * q * Sigma2_tot(mu);
            return ps_nw_q + (1 + k2S2) * std::exp(-k2S2) * ps_w_q;
        }
        double IR_resum_n2(double q, double mu) const {
            double ps_nw_q = ps_nw(q, q_min, q_max);
            double ps_w_q = ps(q, q_min, q_max) - ps_nw_q;
            double k2S2 = q * q * Sigma2_tot(mu);
            return ps_nw_q +
                (1 + k2S2 + 0.5 * k2S2 * k2S2) * std::exp(-k2S2) * ps_w_q;
        }
    public:
        InputPowerSpectrum(
                const std::string& input_ps_filename,
                double input_ps_rescale,
                double q_min,
                double q_max,
                bool ir_resum,
                const IRresumSettings& ir_settings,
                bool rsd,
                double rsd_growth_f
                );
        InputPowerSpectrum(const InputPowerSpectrum& other) = delete;

        double operator()(double q, double mu) const { return evaluate(q, mu); }
};


struct IntegrationInput {
    const double q_min = 0;
    const double q_max = 0;

    const InputPowerSpectrum& ps;

    bool single_hard_limit;

    Vec1D<IntegrandTables> tables_vec;

    /* For power spectrum */
    Vec1D<PowerSpectrumDiagram> ps_diagrams;
    Vec1D<Pair<int>> pair_correlations;
    /* For bispectrum */
    Vec1D<BiSpectrumDiagram> bs_diagrams;
    Vec1D<Triple<int>> triple_correlations;
};


int integrand(
        const int *ndim,
        const double xx[],
        const int *ncomp,
        double ff[],
        void *userdata,
        const int *nvec,
        const int *core
        );


#endif /* ifndef INTEGRAND_HPP */
