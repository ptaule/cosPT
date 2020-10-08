/*
   kernel_evolution.hpp

   Created by Petter Taule on 02.10.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef KERNEL_EVOLUTION_HPP
#define KERNEL_EVOLUTION_HPP

#include "utilities.hpp"
#include "parameters.hpp"
#include "tables.hpp"

int kernel_evolution(
        const int arguments[],
        int kernel_index,
        int n,
        IntegrandTables& tables
        );


void compute_F1(
        double k,
        const EvolutionParameters& params,
        const EtaGrid& eta,
        Vec1D<double> F1_eta_i, /* out */
        Vec1D<double> F1_eta_f  /* out */
        );


#endif /* ifndef KERNEL_EVOLUTION_HPP */
