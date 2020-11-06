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
#include "diagrams.hpp"
#include "tables.hpp"

class IntegrationInput {
    public:
        const double q_min = 0;
        const double q_max = 0;

        const Vec1D<PowerSpectrumDiagram>* ps_diagrams;
        const Vec1D<BiSpectrumDiagram>* bs_diagrams;
        const Vec1D<Pair<int>>* pair_correlations;
        const Vec1D<Triple<int>>* triple_correlations;

        const Interpolation1D& input_ps;
        Vec1D<IntegrandTables>& tables_vec;

        /* For power specturm */
        IntegrationInput(
                double q_min,
                double q_max,
                const Vec1D<PowerSpectrumDiagram>* ps_diagrams,
                const Vec1D<Pair<int>>* pair_correlations,
                const Interpolation1D& input_ps,
                Vec1D<IntegrandTables>& tables_vec
                ) :
            q_min(q_min), q_max(q_max), ps_diagrams(ps_diagrams),
            bs_diagrams(nullptr), pair_correlations(pair_correlations),
            input_ps(input_ps), tables_vec(tables_vec) {}

        /* For bispecturm */
        IntegrationInput(
                double q_min,
                double q_max,
                const Vec1D<BiSpectrumDiagram>* bs_diagrams,
                const Vec1D<Triple<int>>* triple_correlations,
                const Interpolation1D& input_ps,
                Vec1D<IntegrandTables>& tables_vec
                ) :
            q_min(q_min), q_max(q_max), ps_diagrams(nullptr),
            bs_diagrams(bs_diagrams), triple_correlations(triple_correlations),
            input_ps(input_ps), tables_vec(tables_vec) {}
};


namespace ps {
int integrand(
        __attribute__((unused)) const int *ndim,
        const double xx[],
        __attribute__((unused)) const int *ncomp,
        double ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        );
}

namespace bs {
int integrand(
        __attribute__((unused)) const int *ndim,
        const double xx[],
        __attribute__((unused)) const int *ncomp,
        double ff[],
        void *userdata,
        __attribute__((unused)) const int *nvec,
        const int *core
        );
}


class Results {
    private:
        Spectrum spectrum;
        Vec1D<Pair<int>> pair_correlations_;
        Vec1D<Triple<int>> triple_correlations_;
    public:
        Vec1D<double> lin_ps;
        Vec1D<double> non_lin_ps;
        Vec1D<double> errors;

        const Vec1D<Pair<int>>& pair_correlations() const {
            return pair_correlations_;
        }
        const Vec1D<Triple<int>>& triple_correlations() const {
            return triple_correlations_;
        }

        Results(Spectrum spectrum,
                const Vec1D<Pair<int>>& pair_correlations
                ) : spectrum(spectrum), pair_correlations_(pair_correlations)
        {
            lin_ps.resize(pair_correlations_.size());
            non_lin_ps.resize(pair_correlations_.size());
            errors.resize(pair_correlations_.size());
        }

        Results(Spectrum spectrum,
                const Vec1D<Triple<int>>& triple_correlations
                ) : spectrum(spectrum), triple_correlations_(triple_correlations)
        {
            lin_ps.resize(triple_correlations_.size());
            non_lin_ps.resize(triple_correlations_.size());
            errors.resize(triple_correlations_.size());
        }
};


#endif /* ifndef INTEGRAND_HPP */
