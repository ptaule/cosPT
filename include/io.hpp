/*
   io.hpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>

#include "utilities.hpp"
#include "integrand.hpp"


void read_delimited_file(
        const std::string& filename,
        Vec2D<double>& columns
        );


void write_results(
        const Config& cfg,
        const Vec1D<double>& tree_level_result,
        const Vec1D<double>& loop_result,
        const Vec1D<double>& errors
        );

#endif /* ifndef IO_HPP */
