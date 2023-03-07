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

#include "diagrams.hpp"
#include "interpolation.hpp"
#include "utilities.hpp"

class IntegrandTables;
class IRresumSettings;

class InputPowerSpectrum {
    private:
        int loop_order  = 0; /* Which loop order */
        int pt_order = 0; /* At which PT order are we working? To avoid
                              overcounting of soft loops */

        Interpolation1D ps;
        Interpolation1D ps_nw;

        double q_min = 0;
        double q_max = 0;

        bool ir_resum_ = false;
        bool rsd_ = false;
        double rsd_growth_f_ = 0;
        double Sigma2 = 0;
        double delta_Sigma2 = 0;

        double Sigma2_total_rsd(double mu) const {
            double mu2 = mu*mu;
            return (1 + rsd_growth_f_ * mu2 * (2 + rsd_growth_f_)) * Sigma2 +
                rsd_growth_f_ * rsd_growth_f_ * mu2 * (mu2 - 1) * delta_Sigma2;
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

        std::function<double(double)> Sigma2_tot;
        std::function<double(double, double)> evaluate;
    public:
        InputPowerSpectrum(
                const std::string& input_ps_filename,
                double input_ps_rescale,
                double q_min,
                double q_max,
                bool ir_resum,
                const IRresumSettings& ir_settings,
                int loop_order,
                int pt_order,
                bool rsd,
                double rsd_growth_f
                );
        InputPowerSpectrum(const InputPowerSpectrum& other) = delete;

        /* Power spectrum that is input to loop calculations */
        double operator()(double q, double mu) const { return evaluate(q, mu); }

        /* Tree-level power spectrum (depends also on pt_order) */
        double tree_level(double q, double mu) const;

        /* Getters */
        bool ir_resum() const { return ir_resum_; }
        bool rsd() const { return rsd_; }
        double rsd_growth_f() const { return rsd_growth_f_; }
};


struct IntegrationInput {
    const InputPowerSpectrum& ps;

    const double q_min = 0;
    const double q_max = 0;

    bool single_hard_limit;

    Vec1D<IntegrandTables> tables_vec;

    /* For power spectrum */
    Vec1D<PowerSpectrumDiagram> ps_diagrams;
    Vec1D<Pair<int>> pair_correlations;
    /* For bispectrum */
    Vec1D<BiSpectrumDiagram> bs_diagrams;
    Vec1D<Triple<int>> triple_correlations;

    IntegrationInput(
            const InputPowerSpectrum& ps,
            double q_min,
            double q_max,
            bool single_hard_limit
            ) : ps(ps), q_min(q_min), q_max(q_max),
    single_hard_limit(single_hard_limit) {}
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
