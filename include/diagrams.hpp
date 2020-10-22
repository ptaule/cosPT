/*
   diagrams.hpp

   Created by Petter Taule on 30.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
   */

#ifndef DIAGRAMS_HPP
#define DIAGRAMS_HPP

#include <vector>
#include <iosfwd>

#include "utilities.hpp"
#include "parameters.hpp"

class ArgumentConfiguration {
    public:
        int kernel_index;
        Vec1D<int> args;
};

class PowerSpectrumDiagram {
    private:
        const LoopParameters& loop_params;

        void compute_rearrangements(int n_loops);
        void compute_sign_flips();

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(int n_coeffs, int a, int b);
    public:
        int m, l, r;

        int diagram_factor;     /* Topological multiplicative diagram factor */

        /* Table of rearrangements */
        Vec2D<int> rearrangements;
        /* Table of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Argument configuration for each rearrangement and sign setup
         * (n_rearrangements x n_sign_configs possibilities) */
        Vec2D<ArgumentConfiguration> arg_configs_l;
        Vec2D<ArgumentConfiguration> arg_configs_r;

        PowerSpectrumDiagram(
                const LoopParameters& params,
                int m, int l, int r
                );

        void print_diagram_tags(std::ostream& out) const;
        void print_argument_configuration(
                std::ostream& out,
                int a,
                int b
                ) const;
};


class BiSpectrumDiagram {
    private:
        const LoopParameters& loop_params;

        /* Is diagram closed (n_ab, n_bc, n_ca > 0) */
        bool overall_loop;

        /* Number of sign_flips for connecting lines */
        int n_connecting_loops_ab;
        int n_connecting_loops_bc;
        int n_connecting_loops_ca;

        void compute_rearrangements(int n_loops);
        Vec2D<bool> connecting_line_sign_flips(int n_connecting_loops) const;
        Vec2D<bool> compute_sign_flips(
                const Vec2D<bool>& sign_configs_ab,
                const Vec2D<bool>& sign_configs_bc,
                const Vec2D<bool>& sign_configs_ca
                );

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(
              int n_coeffs,
              int rearr_idx,
              int sign_idx,
              int overall_loop_idx
              );

    public:
        int n_ab, n_bc, n_ca;
        int n_a, n_b, n_c;

        int diagram_factor;     /* Topological multiplicative diagram factor */

        /* Table of rearrangements */
        Vec2D<int> rearrangements;
        /* Tables of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Argument configuration for each rearrangement, sign setup and
         * (potential) overall loop assosiation
         * (n_rearrangements x n_sign_configs x 3 possibilities) */
        Vec3D<ArgumentConfiguration> arg_configs_a;
        Vec3D<ArgumentConfiguration> arg_configs_b;
        Vec3D<ArgumentConfiguration> arg_configs_c;

        BiSpectrumDiagram(
                const LoopParameters& loop_params,
                int n_ab, int n_bc, int n_ca,
                int n_a, int n_b, int n_c
                );

        bool has_overall_loop() const {return overall_loop;}

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
