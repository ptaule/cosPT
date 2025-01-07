#ifndef INTEGRAND_HPP
#define INTEGRAND_HPP

#include "diagrams.hpp"
#include "utilities.hpp"

class IntegrandTables;
class InputPowerSpectrum;


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
