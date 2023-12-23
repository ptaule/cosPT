/*
   tree_level.hpp

   Created by Petter Taule on 19.12.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef TREE_LEVEL_HPP
#define TREE_LEVEL_HPP

#include "utilities.hpp"

class Interpolation1D;
class IntegrandTables;
class InputPowerSpectrum;
class EvolutionParameters;
class EtaGrid;


namespace ps {
void tree_level(
    double k_a,
    Dynamics dynamics,
    const InputPowerSpectrum& ps,
    const EtaGrid& eta_grid,
    const EvolutionParameters& ev_params,
    const Vec1D<Pair<int>>& pair_correlations,
    Vec1D<double>& results /* out */
);
} /* namespace ps */


namespace bs {
void tree_level(
        IntegrandTables& tables,
        const InputPowerSpectrum& ps,
        const Vec1D<Triple<int>>& triple_correlations,
        Vec1D<double>& results /* out */
        );
} /* namespace bs */


#endif /* ifndef TREE_LEVEL_HPP */
