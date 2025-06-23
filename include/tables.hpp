#ifndef TABLES_HPP
#define TABLES_HPP

#include <array>
#include <stdexcept>

#include "utilities.hpp"

class LoopParameters;
class EvolutionParameters;
class OmegaEigenspace;

struct IntegrationVariables {
    Vec1D<double> magnitudes; /* Loop momenta magnitudes */
    Vec1D<double> cos_theta;  /* Cosine of polar angles of the loop momenta */
    Vec1D<double> phi;        /* Azimutal angles */
    double mu_los;            /* RSD: Cosine of angle between k and L.o.S. */

    IntegrationVariables(std::size_t n_loops) : mu_los(0) {
        magnitudes.assign(n_loops,0);
        cos_theta.assign(n_loops,0);
        phi.assign(n_loops,0);
    }
};


struct SPTKernel {
    std::array<double, EDS_SPT_COMPONENTS> values;
    bool computed = false;
};


struct Kernel {
    Vec2D<double> values;
    bool computed = false;
};


struct RSDKernel {
    double value;
    bool computed = false;
};


class SumTable {
    private:
        const int zero_label;
        const std::size_t n_coeffs;

        Vec2D<int> sum_table;

        void check_result(int res) const;
        int convert_and_sum(int a, int b);
    public:
        SumTable(const LoopParameters& loop_params);

        int operator()(int a, int b) const {
            int result = sum_table.at(static_cast<size_t>(a))
                                  .at(static_cast<size_t>(b));
#if DEBUG > 0
            check_result(result);
#endif
            return result;
        }

        int operator()(const int labels[], std::size_t size) const;
};


class EtaGrid {
    private:
        double eta_ini_             = 0;
        double eta_fin_             = 0;
        std::size_t time_steps_     = 0;
        std::size_t pre_time_steps_ = 0;
        double eta_asymp_           = 0;

        Vec1D<double> grid_;
    public:
        EtaGrid() = default;
        EtaGrid(
                double eta_ini,
                double eta_fin,
                std::size_t time_steps,
                std::size_t pre_time_steps,
                double eta_asymp
                );

        double eta_ini() const {return eta_ini_;}
        double eta_fin() const {return eta_fin_;}
        std::size_t time_steps() const {return time_steps_;}
        std::size_t pre_time_steps() const {return pre_time_steps_;}
        double eta_asymp() const {return eta_asymp_;}

        const Vec1D<double>& grid() const {return grid_;}

        const double& operator[](std::size_t i) const {return grid_[i];}
        const double& at(size_t i) const {return grid_.at(i);}
};


std::ostream& operator<<(std::ostream& out, const EtaGrid& eta_grid);


/* Tables inside integration, one instance for each integration thread */
class IntegrandTables {
    private:
        double k_a    = 0.0;
        double k_b    = 0.0; // For bispectrum
        double cos_ab = 0.0; // For bispectrum

        double rsd_f; /* Growth factor (at observation redshift) for L.o.S. */

        /* Helper vectors for compute_scalar_products() */
        Vec1D<int> a_coeffs;
        Vec1D<int> b_coeffs;

        Vec2D<double> bare_dot_prod;
        /* N_COEFFS x N_COEFFS, dot products between external/loop wavenumbers */
        Vec2D<double> composite_dot_prod;
        /* N_CONFIGS x N_CONFIGS, dot products between composite wavenumbers */
        Vec2D<double> alpha_;                /* N_CONFIGS x N_CONFIGS */
        Vec2D<double> beta_;                 /* N_CONFIGS x N_CONFIGS */

        Vec1D<double> bare_los_projection_;
        /* N_COEFFS dot products between external/loop wavenumber and L.o.S */
        Vec1D<double> composite_los_projection_;
        /* N_CONFIGS dot products between composite wavenumber and L.o.S */

        void reset_spt_kernels();
        void reset_kernels();
        void reset_rsd_kernels();

        void compute_bare_los_proj();
        void compute_composite_los_proj();

        void compute_bare_dot_prod();

        void compute_composite_dot_prod();
        void compute_alpha_beta();
    public:
        const LoopParameters& loop_params;
        const SumTable& sum_table;

        const EvolutionParameters& ev_params;
        const EtaGrid& eta_grid;
        const OmegaEigenspace& omega_eigenspace;

        IntegrationVariables vars;

        Vec1D<SPTKernel> spt_kernels;
        Vec1D<Kernel> kernels;

        Vec1D<RSDKernel> rsd_kernels;
        Vec2D<RSDKernel> vel_power_kernels;

        IntegrandTables() = delete;
        IntegrandTables(
                double k_a,
                double k_b,
                double cos_ab,
                double rsd_growth_f,
                const LoopParameters& loop_params,
                const SumTable& sum_table,
                const EvolutionParameters& ev_params,
                const EtaGrid& eta_grid,
                const OmegaEigenspace& omega_eigenspace
                );

        const Vec2D<double>& bare_dot_products() const { return bare_dot_prod; }
        const Vec2D<double>& composite_dot_products() const { return composite_dot_prod; }
        const Vec2D<double>& alpha() const {return alpha_;}
        const Vec2D<double>& beta() const {return beta_;}

        const Vec1D<double>& bare_los_projection() const { return bare_los_projection_; }
        const Vec1D<double>& composite_los_projection() const { return composite_los_projection_; }

        double get_k_a() const { return k_a; }
        double rsd_growth_f() const { return rsd_f; }

        void reset();
        void compute_tables();
};


#endif /* ifndef TABLES_HPP */
