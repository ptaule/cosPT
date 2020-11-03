/*
   diagrams.hpp

   Created by Petter Taule on 30.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
   */

#ifndef DIAGRAMS_HPP
#define DIAGRAMS_HPP

#include <vector>
#include <iosfwd>
#include <cmath>

#include "utilities.hpp"
#include "parameters.hpp"

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
        void kernel_arguments(int rearr_idx, int sign_idx);
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
                int rearr_idx,
                int sign_idx,
                const Vec2D<double>& scalar_products
                ) const
        {
            int q_m1_label = q_m1_labels.at(rearr_idx).at(sign_idx);
            return std::sqrt(scalar_products.at(q_m1_label).at(q_m1_label));
        }

        /* Heaviside theta for q_m1 and first connecting loop */
        double heaviside_theta(
                double q_m1,
                int rearr_idx,
                const Vec1D<double>& loop_magnitudes
                ) const
        {
            const Vec1D<int>& rearr = rearrangements.at(rearr_idx);
            if (m == 1) return 1;
            // Heaviside-theta (q_m1 - Q1(rearr))
            if (q_m1 <= loop_magnitudes.at(rearr.at(0))) return 0;
#if DEBUG >= 1
            /* Check that the heaviside-theta (Q2(rearr) - Q3(rearr)) etc. are
             * satisfied by (reparametrized) momenta from CUBA */
            for (int i = 3; i <= m; ++i) {
                if (loop_magnitudes.at(rearr.at(i - 3)) <
                    loop_magnitudes.at(rearr.at(i - 2))
                    ) {
                    throw(std::logic_error(
                        "Heaviside theta: Q" + std::to_string(rearr.at(i - 3) + 1) +
                        " < Q" + std::to_string(rearr.at(i - 2) + 1) + "."));
                }
            }
#endif
            return m;
        }

        const ArgumentConfiguration& get_arg_config_l(int rearr_idx, int
                sign_idx) const {
            return arg_configs_l.at(rearr_idx).at(sign_idx);
        }
        const ArgumentConfiguration& get_arg_config_r(int rearr_idx, int
                sign_idx) const {
            return arg_configs_r.at(rearr_idx).at(sign_idx);
        }
#undef at

        void print_diagram_tags(std::ostream& out) const;
        void print_argument_configuration(
                std::ostream& out,
                int a,
                int b
                ) const;
};


class BiSpectrumDiagram {
    private:
        /* Is diagram closed (n_ab, n_bc, n_ca > 0) */
        bool overall_loop_;
        /* Topological multiplicative diagram factor */
        int diagram_factor_;

        /* Number of connecting loops */
        int n_connecting_loops_ab;
        int n_connecting_loops_bc;
        int n_connecting_loops_ca;

        const LoopParameters& loop_params;

        /* Table of rearrangements */
        Vec2D<int> rearrangements;
        /* Tables of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Matrices of q_ab1, q_bc1 and q_ca1 labels, for each combination of
         * rearr_idx, sign_idx and overall_loop_idx */
        Vec3D<int> q_ab1_labels;
        Vec3D<int> q_bc1_labels;
        Vec3D<int> q_ca1_labels;

        /* Argument configuration for each rearrangement, sign setup and
         * (potential) overall loop assosiation
         * (n_rearrangements x n_sign_configs x 3 possibilities) */
        Vec3D<ArgumentConfiguration> arg_configs_a;
        Vec3D<ArgumentConfiguration> arg_configs_b;
        Vec3D<ArgumentConfiguration> arg_configs_c;

        void compute_rearrangements(int n_loops);
        Vec2D<bool> connecting_line_sign_flips(int n_connecting_loops) const;
        Vec2D<bool> compute_sign_flips(
                const Vec2D<bool>& sign_configs_ab,
                const Vec2D<bool>& sign_configs_bc,
                const Vec2D<bool>& sign_configs_ca
                );

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(int rearr_idx, int sign_idx, int overall_loop_idx);

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

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif
        /* q_ab1 is the momentum of the connecting line ab whose loop momentum
         * is evaluated by the momentum conserving delta function, similarly
         * for q_bc1 and q_ca1 */
        double q_ab1(
                int rearr_idx,
                int sign_idx,
                int overall_loop_idx,
                const Vec2D<double>& scalar_products
                ) const
        {
            int q_ab1_label =
                q_ab1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
            return std::sqrt(scalar_products.at(q_ab1_label).at(q_ab1_label));
        }

        double q_bc1(
                int rearr_idx,
                int sign_idx,
                int overall_loop_idx,
                const Vec2D<double>& scalar_products
                ) const
        {
            int q_bc1_label =
                q_bc1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
            return std::sqrt(scalar_products.at(q_bc1_label).at(q_bc1_label));
        }

        double q_ca1(
                int rearr_idx,
                int sign_idx,
                int overall_loop_idx,
                const Vec2D<double>& scalar_products
                ) const
        {
            int q_ca1_label =
                q_ca1_labels.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
            return std::sqrt(scalar_products.at(q_ca1_label).at(q_ca1_label));
        }

        const ArgumentConfiguration& get_arg_config_a(int rearr_idx, int
                sign_idx, int overall_loop_idx) const {
            return arg_configs_a.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
        }
        const ArgumentConfiguration& get_arg_config_b(int rearr_idx, int
                sign_idx, int overall_loop_idx) const {
            return arg_configs_b.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
        }
        const ArgumentConfiguration& get_arg_config_c(int rearr_idx, int
                sign_idx, int overall_loop_idx) const {
            return arg_configs_c.at(rearr_idx).at(sign_idx).at(overall_loop_idx);
        }
#undef at

        void print_diagram_tags(std::ostream& out) const;
        void print_argument_configuration(
                std::ostream& out,
                int rearr_idx,
                int sign_idx,
                int overall_loop_idx = 0
                ) const;
};

std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram);
std::ostream& operator<<(std::ostream& out, const BiSpectrumDiagram& diagram);

namespace ps {
    Vec1D<PowerSpectrumDiagram> construct_diagrams(const LoopParameters&
            loop_params);
}

namespace bs {
    Vec1D<BiSpectrumDiagram> construct_diagrams(const LoopParameters&
            loop_params);
}

#endif /* ifndef DIAGRAMS_HPP */
