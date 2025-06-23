#ifndef DIAGRAMS_HPP
#define DIAGRAMS_HPP

#include <cmath>
#include <iosfwd>
#include <string>

#include "utilities.hpp"

class LoopParameters;

class ArgumentConfiguration {
    public:
        int kernel_index;
        Vec1D<int> args;
};

class PowerSpectrumDiagram {
    private:
        int diagram_factor_; /* Topological multiplicative diagram factor */

        const LoopParameters& loop_params;

        /* Table of rearrangements */
        Vec2D<int> rearrangements;
        /* Table of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Matrix of q_m1 labels, for each combination of rearr_idx and sign_idx */
        Vec2D<int> q_m1_labels;

        /* Argument configuration for each rearrangement and sign setup
         * (n_rearrangements x n_sign_configs possibilities) */
        Vec2D<ArgumentConfiguration> arg_configs_l;
        Vec2D<ArgumentConfiguration> arg_configs_r;

        void compute_rearrangements(int n_loops);
        void compute_sign_flips();

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(std::size_t rearr_idx, std::size_t sign_idx);
    public:
        const int m, l, r;

        PowerSpectrumDiagram(
                const LoopParameters& params,
                int m, int l, int r
                );

        int diagram_factor() const {return diagram_factor_;}
        size_t n_rearrangements() const {return rearrangements.size();}
        size_t n_sign_configs() const {return sign_configs.size();}

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif
        /* q_m1 is the momentum of the connecting line whose loop momentum is
         * evaluated by the momentum conserving delta function */
        double q_m1(
                std::size_t rearr_idx,
                std::size_t sign_idx,
                const Strided2DVec<double>& dot_products
                ) const
        {
            std::size_t q_m1_label =
                static_cast<std::size_t>(q_m1_labels.at(rearr_idx).at(sign_idx));
            return std::sqrt(dot_products(q_m1_label, q_m1_label));
        }

        /* Heaviside theta for q_m1 and first connecting loop */
        int heaviside_theta(
                double q_m1,
                std::size_t rearr_idx,
                const Vec1D<double>& loop_magnitudes
                ) const
        {
            const Vec1D<int>& rearr = rearrangements.at(rearr_idx);
            if (m == 1) return 1;
            // Heaviside-theta (q_m1 - Q1(rearr))
            if (q_m1 <= loop_magnitudes.at(static_cast<size_t>(rearr.at(0)))) {
                return 0;
            }
#if DEBUG >= 1
            /* Check that the heaviside-theta (Q2(rearr) - Q3(rearr)) etc. are
             * satisfied by (reparametrized) momenta from CUBA */
            for (size_t i = 3; i <= static_cast<size_t>(m); ++i) {
                if (loop_magnitudes.at(static_cast<size_t>(rearr.at(i - 3))) <
                    loop_magnitudes.at(static_cast<size_t>(rearr.at(i - 2)))
                    ) {
                    throw(std::logic_error(
                        "Heaviside theta: Q" + std::to_string(rearr.at(i - 3) + 1) +
                        " < Q" + std::to_string(rearr.at(i - 2) + 1) + "."));
                }
            }
#endif
            return m;
        }

        const ArgumentConfiguration & get_arg_config_l(
                std::size_t rearr_idx,
                std::size_t sign_idx
                ) const {
            return arg_configs_l.at(rearr_idx).at(sign_idx);
        }
        const ArgumentConfiguration& get_arg_config_r(
                std::size_t rearr_idx,
                std::size_t sign_idx
                ) const {
            return arg_configs_r.at(rearr_idx).at(sign_idx);
        }
#undef at

        std::string tags() const {
            return ("(m,l,r) = ("
                + std::to_string(m) + ","
                + std::to_string(l) + ","
                + std::to_string(r) + ")");
        }
        std::string argument_configuration(std::size_t a,std::size_t b) const;
};


class BiSpectrumDiagram {
    private:
        /* Is diagram closed (n_ab, n_bc, n_ca > 0) */
        bool overall_loop_;
        /* Topological multiplicative diagram factor */
        int diagram_factor_;

        /* Number of connecting loops */
        Triple<int> n_connecting_loops;

        const LoopParameters& loop_params;

        /* Table of rearrangements */
        Vec2D<int> rearrangements;
        /* Tables of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;
        /* Matrices of q_xy1 = (q_ab1, q_bc1, q_ca1) labels, for each combination of
         * rearr_idx, sign_idx and overall_loop_idx */
        Vec3D<Triple<int>> q_xy1_labels;
        /* Vectors of q_xy = (q_ab, q_bc, q_ca) (sum of connecting lines)
         * labels, for each rearr_idx and overall_loop_idx */
        Vec2D<Triple<int>> q_xy_labels;

        /* Argument configuration for each rearrangement, sign setup and
         * (potential) overall loop assosiation
         * (n_rearrangements x n_sign_configs x 3 possibilities) */
        Vec3D<Triple<ArgumentConfiguration>> arg_configs;

        void compute_rearrangements(int n_loops);
        Vec2D<bool> connecting_line_sign_flips(int connecting_loops) const;
        void compute_sign_flips(const Triple<Vec2D<bool>>& sign_configs_xy);

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(
                std::size_t rearr_idx,
                std::size_t sign_idx,
                std::size_t overall_loop_idx
                );

    public:
        const int n_ab, n_bc, n_ca;
        const int n_a, n_b, n_c;

        BiSpectrumDiagram(
                const LoopParameters& loop_params,
                int n_ab, int n_bc, int n_ca,
                int n_a, int n_b, int n_c
                );

        bool overall_loop() const {return overall_loop_;}
        int diagram_factor() const {return diagram_factor_;}
        size_t n_rearrangements() const {return rearrangements.size();}
        size_t n_sign_configs() const {return sign_configs.size();}

        void connecting_lines_factors(
                std::size_t rearr_idx,
                std::size_t sign_idx,
                std::size_t overall_loop_idx,
                const Vec1D<double>& loop_magnitudes,
                const Strided2DVec<double>& dot_products,
                Triple<double>& q_xy1, /* out */
                int& heaviside_theta   /* out */
                ) const;

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif
        /* Get the magnitude of the first (rearranged) loop momentum */
        double q1_magnitude(
                std::size_t rearr_idx,
                const Vec1D<double>& loop_magnitudes
                ) const
        {
            return loop_magnitudes.at(
                static_cast<size_t>(rearrangements.at(rearr_idx).at(0)));
        }

        const Triple<ArgumentConfiguration>& get_arg_config(
                std::size_t rearr_idx,
                std::size_t sign_idx,
                std::size_t overall_loop_idx
                ) const {
            return arg_configs.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
        }
#undef at

        std::string tags() const;
        std::string argument_configuration(
                std::size_t rearr_idx,
                std::size_t sign_idx,
                std::size_t overall_loop_idx = 0
                ) const;
};

std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram);
std::ostream& operator<<(std::ostream& out, const BiSpectrumDiagram& diagram);

namespace ps {
Vec1D<PowerSpectrumDiagram> construct_diagrams(const LoopParameters& loop_params);
}

namespace bs {
Vec1D<BiSpectrumDiagram> construct_diagrams(const LoopParameters& loop_params);
}

#endif /* ifndef DIAGRAMS_HPP */
