/*
   bispectrum_tree_level.hpp

   Created by Petter Taule on 19.12.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef BISPECTRUM_TREE_LEVEL_HPP
#define BISPECTRUM_TREE_LEVEL_HPP

#include "utilities.hpp"

class Interpolation1D;
class IntegrandTables;

void tree_level_bispectrum(
        IntegrandTables& tables,
        const Interpolation1D& input_ps,
        const Vec1D<Triple<int>>& triple_correlations,
        Vec1D<double>& results /* out */
        );

#endif /* ifndef BISPECTRUM_TREE_LEVEL_HPP */
