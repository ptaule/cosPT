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
#include "tables.hpp"

class ArgumentConfiguration {
    public:
        short int kernel_index;
        Vec1D<short int> args;
};

class PowerSpectrumDiagram {
    private:
        const Settings& settings;

        void compute_rearrangements(short int n_loops);
        void compute_sign_flips();

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(short int n_coeffs, short int a, short int b);
    public:
        short int m;
        short int l;
        short int r;

        short int diagram_factor;     /* Topological multiplicative diagram factor */

        short int n_rearrangements;   /* Number of rearrangements of loop momenta  */
        short int n_sign_configs;     /* Number of sign flips                      */
        /* Table of rearrangements */
        Vec2D<short int> rearrangements;
        /* Table of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Argument configuration for each rearrangement and sign setup
         * (n_rearrangements x n_sign_configs possibilities) */
        Vec2D<ArgumentConfiguration> arg_configs_l;
        Vec2D<ArgumentConfiguration> arg_configs_r;

        PowerSpectrumDiagram(
                const Settings& settings,
                short int m,
                short int l,
                short int r
                );

        void print_diagram_tags(std::ostream& out) const;
        void print_argument_configuration(
                std::ostream& out,
                short int a,
                short int b
                ) const;
};


class BiSpectrumDiagram {
    private:
        const Settings& settings;

        /* Is diagram closed (n_ab, n_bc, n_ca > 0) */
        bool overall_loop;

        /* Number of sign_flips for connecting lines */
        short int n_connecting_loops_ab;
        short int n_connecting_loops_bc;
        short int n_connecting_loops_ca;

        void compute_rearrangements(short int n_loops);
        Vec2D<bool> connecting_line_sign_flips(short int n_connecting_loops) const;
        Vec2D<bool> compute_sign_flips(
                const Vec2D<bool>& sign_configs_ab,
                const Vec2D<bool>& sign_configs_bc,
                const Vec2D<bool>& sign_configs_ca
                );

        /* Computes argument configurations if rearrangement and sign
         * configurations are set */
        void kernel_arguments(
              short int n_coeffs,
              short int rearr_idx,
              short int sign_idx,
              short int overall_loop_idx
              );

    public:
        short int n_ab, n_bc, n_ca;
        short int n_a, n_b, n_c;

        short int diagram_factor;     /* Topological multiplicative diagram factor */

        short int n_rearrangements;   /* Number of rearrangements of loop momenta  */

        /* Table of rearrangements */
        Vec2D<short int> rearrangements;
        /* Tables of sign flips, true <-> +1, false <-> -1 */
        Vec2D<bool> sign_configs;

        /* Argument configuration for each rearrangement, sign setup and
         * (potential) overall loop assosiation
         * (n_rearrangements x n_sign_configs x 3 possibilities) */
        Vec3D<ArgumentConfiguration> arg_configs_a;
        Vec3D<ArgumentConfiguration> arg_configs_b;
        Vec3D<ArgumentConfiguration> arg_configs_c;

        BiSpectrumDiagram(
                const Settings& settings,
                short int n_ab, short int n_bc, short int n_ca,
                short int n_a, short int n_b, short int n_c
                );

        void print_diagram_tags(std::ostream& out) const;
        /* void print_argument_configuration( */
        /*         std::ostream& out, */
        /*         short int a, */
        /*         short int b */
        /*         ) const; */
};

std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram);
std::ostream& operator<<(std::ostream& out, const BiSpectrumDiagram& diagram);

Vec1D<PowerSpectrumDiagram> construct_ps_diagrams(const Settings& settings);
Vec1D<BiSpectrumDiagram>    construct_bs_diagrams(const Settings& settings);

#endif /* ifndef DIAGRAMS_HPP */
