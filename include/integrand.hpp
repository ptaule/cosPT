/*
   integrand.hpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef INTEGRAND_HPP
#define INTEGRAND_HPP

#include <utility>

#include "utilities.hpp"
#include "interpolation.hpp"

class IntegrandTables;
class PowerSpectrumDiagram;
class BiSpectrumDiagram;

struct IntegrationInput {
    const double q_min = 0;
    const double q_max = 0;

    bool single_hard_limit;

    Interpolation1D input_ps;
    Vec1D<IntegrandTables> tables_vec;

    /* For power spectrum */
    Vec1D<PowerSpectrumDiagram> ps_diagrams;
    Vec1D<Pair<int>> pair_correlations;
    /* For bispectrum */
    Vec1D<BiSpectrumDiagram> bs_diagrams;
    Vec1D<Triple<int>> triple_correlations;

    IntegrationInput(double q_min, double q_max, bool single_hard_limit)
    : q_min(q_min), q_max(q_max), single_hard_limit(single_hard_limit) {}
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
