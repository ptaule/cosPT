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
        void kernel_arguments(short int n_loops, short int a, short int b);
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

std::ostream& operator<<(std::ostream& out, const PowerSpectrumDiagram& diagram);

Vec1D<PowerSpectrumDiagram> construct_diagrams(const Settings& settings);

#endif /* ifndef DIAGRAMS_HPP */
