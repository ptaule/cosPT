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

using Correlation = std::pair<int, int>;

class IntegrationInput {
    public:
        const double q_min = 0;
        const double q_max = 0;

        const Vec1D<PowerSpectrumDiagram>* ps_diagrams;
        const Vec1D<BiSpectrumDiagram>* bs_diagrams;
        const Vec1D<Correlation>& correlations;
        const Interpolation1D& input_ps;

        Vec1D<IntegrandTables>& tables_vec;

        IntegrationInput(
                double q_min,
                double q_max,
                const Vec1D<PowerSpectrumDiagram>* ps_diagrams,
                const Vec1D<Correlation>& correlations,
                const Interpolation1D& input_ps,
                Vec1D<IntegrandTables>& tables_vec
                ) :
            q_min(q_min), q_max(q_max), ps_diagrams(ps_diagrams),
            bs_diagrams(nullptr), correlations(correlations),
            input_ps(input_ps), tables_vec(tables_vec) {}

        IntegrationInput(
                double q_min,
                double q_max,
                const Vec1D<BiSpectrumDiagram>* bs_diagrams,
                const Vec1D<Correlation>& correlations,
                const Interpolation1D& input_ps,
                Vec1D<IntegrandTables>& tables_vec
                ) :
            q_min(q_min), q_max(q_max), ps_diagrams(nullptr),
            bs_diagrams(bs_diagrams), correlations(correlations),
            input_ps(input_ps), tables_vec(tables_vec) {}
};


namespace ps {
    void integrand(
            const IntegrationInput& input,
            IntegrandTables& tables,
            Vec1D<double>& results
            );
}

namespace bs {
    void integrand(
            const IntegrationInput& input,
            IntegrandTables& tables,
            Vec1D<double>& results
            );
}


class Results {
    private:
        const Vec1D<Correlation>& correlations;
    public:
        Vec1D<double> lin_ps;
        Vec1D<double> non_lin_ps;
        Vec1D<double> errors;

        const Vec1D<Correlation>& get_correlations() const {
            return correlations;
        }

        Results(const Vec1D<Correlation>& correlations)
            : correlations(correlations)
        {
            lin_ps.resize(correlations.size());
            non_lin_ps.resize(correlations.size());
            errors.resize(correlations.size());
        }
};


std::ostream& operator<<(std::ostream&, const Correlation& correlation);

#endif /* ifndef INTEGRAND_HPP */
