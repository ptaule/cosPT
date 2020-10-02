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


void read_columns_from_file(
        const std::string& filename, /* in, name of file to be read */
        unsigned int n_columns,      /* in, number of columns       */
        Vec2D<double>& columns       /* out, columns                */
        );


void read_data_grid_from_file(
        const std::string& filename,
        Vec1D<double>& data,
        unsigned int n_rows,
        unsigned int n_columns
        );

void write_results(
        const std::string& output_file,
        const std::string& input_ps_file,
        const std::string& description,
        int n_loops,
        double cuba_epsabs,
        double cuba_epsrel,
        double cuba_maxevals,
        double k,
        double q_min,
        double q_max,
        const Results& results
        );

#endif /* ifndef IO_HPP */
