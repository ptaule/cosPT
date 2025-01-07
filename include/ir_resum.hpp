#ifndef IR_RESUM_HPP
#define IR_RESUM_HPP

#include <cstddef>
#include <functional>
#include <string>

#include <gsl/gsl_integration.h>

#include "interpolation.hpp"

class Interpolation1D;


void DST_II(std::vector<double>& data);
void DST_III(std::vector<double>& data);


struct IRresumSettings {
    /* Wiggly/non-wiggly split */
    double k_max = 10;
    int N_power = 16;
    int N_left = 120;
    int N_right = 240;

    /* Damping factor */
    double k_min = 1e-4;
    double k_s = 0.2;
    double k_osc = 1.0/110.0;
    std::size_t integration_sub_regions = 10000;
    double integration_atol = 0;
    double integration_rtol = 1e-6;
    int integration_key = GSL_INTEG_GAUSS61;

    IRresumSettings(double k_s, double k_osc) : k_s(k_s), k_osc(k_osc) {}
};


class InputPowerSpectrum {
    private:
        int loop_order  = 0; /* Which loop order */
        int pt_order = 0; /* At which PT order are we working? To avoid
                              overcounting of soft loops */

        Interpolation1D ps;
        Interpolation1D ps_nw;

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
            double ps_nw_q = ps_nw(q);
            double ps_w_q = ps(q) - ps_nw_q;
            double k2S2 = q * q * Sigma2_tot(mu);
            return ps_nw_q + std::exp(-k2S2) * ps_w_q;
        }
        double IR_resum_n1(double q, double mu) const {
            double ps_nw_q = ps_nw(q);
            double ps_w_q = ps(q) - ps_nw_q;
            double k2S2 = q * q * Sigma2_tot(mu);
            return ps_nw_q + (1 + k2S2) * std::exp(-k2S2) * ps_w_q;
        }
        double IR_resum_n2(double q, double mu) const {
            double ps_nw_q = ps_nw(q);
            double ps_w_q = ps(q) - ps_nw_q;
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

        /* Integral of loop-level power spectrum from a to b (not exceeding
         * interpolation limits) */
        double integral(
            double a,
            double b,
            std::size_t integration_sub_regions = 10000,
            double integration_atol = 0,
            double integration_rtol = 1e-6,
            int integration_key = GSL_INTEG_GAUSS61
        ) const;

        /* Getters */
        bool ir_resum() const { return ir_resum_; }
        bool rsd() const { return rsd_; }
        double rsd_growth_f() const { return rsd_growth_f_; }
};


/* Remove BAO wiggles. Default parameters taken from CLASS-PT, i.e. 2004.10607 */
void remove_BAO_wiggles(
        const Interpolation1D& ps,
        Interpolation1D& ps_nw,
        const IRresumSettings& settings
        );


void compute_ir_damping(
        const Interpolation1D& ps_nw,
        double& Sigma2,       /* out */
        double& delta_Sigma2, /* out */
        const IRresumSettings& settings
        );


#endif /* ifndef IR_RESUM_HPP */
