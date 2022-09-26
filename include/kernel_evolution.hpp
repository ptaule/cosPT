/*
   kernel_evolution.hpp

   Created by Petter Taule on 02.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef KERNEL_EVOLUTION_HPP
#define KERNEL_EVOLUTION_HPP

#include "utilities.hpp"

class IntegrandTables;
class EvolutionParameters;
class EtaGrid;

int compute_gen_kernels(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


void compute_F1(
        double k,
        const EvolutionParameters& ev_params,
        const EtaGrid& eta_grid,
        Vec1D<double>& F1_eta_ini, /* out */
        Vec1D<double>& F1_eta_fin  /* out */
        );


#endif /* ifndef KERNEL_EVOLUTION_HPP */
