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

void tree_level_bispectrum(
        IntegrandTables& tables,
        const InputPowerSpectrum& ps,
        const Vec1D<Triple<int>>& triple_correlations,
        Vec1D<double>& results /* out */
        );

#endif /* ifndef TREE_LEVEL_HPP */
