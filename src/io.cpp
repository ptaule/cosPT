/*
   io.cpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "../include/utilities.hpp"
#include "../include/integrand.hpp"
#include "../include/version.hpp"
#include "../include/io.hpp"

#define RESERVE_SIZE 200

using std::size_t;
using std::setw;

/* Read file with n_columns number of columns, if there are lines that does not
 * have this number of columns, throw error. Skips empty lines and lines
 * beginning with '#' */
void read_columns_from_file(
        const std::string& filename,  /* in, name of file to be read */
        unsigned int n_columns, /* in, number of columns       */
        Vec2D<double>& columns        /* out, columns                */
        )
{
    std::ifstream input(filename);
    std::string line;

    if (input.fail()){
        throw(std::invalid_argument("Could not find " + filename + "."));
    }

    columns.resize(n_columns);
    for (unsigned int i = 0; i < n_columns; ++i) {
        columns.at(i).reserve(RESERVE_SIZE);
    }

    while (getline(input, line)) {
        /* Ignore empty lines or lines beginning with # */
        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        std::stringstream ss(line);
        double value;

        // Column counter
        unsigned int i = 0;
        while (ss >> value) {
            if (i < n_columns) {
                columns.at(i).push_back(value);
            }
            else {
                throw(std::runtime_error("Found line with number of columns "
                                         "not equal to n_columns in " +
                                         filename));
            }
            i++;
        }
        if (i != n_columns) {
            throw(std::runtime_error(
                "Found line with number of columns not equal to n_columns in " +
                filename));
        }
    }
    input.close();
}


/* Reads data from filename, using indexing scheme:
 * index = row + column * n_rows
 * If the number of rows/columns found does not equal n_rows/n_columns, an
 * error is thrown */
void read_data_grid_from_file(
        const std::string& filename,
        Vec1D<double>& data,
        unsigned int n_rows,
        unsigned int n_columns
        )
{
    std::ifstream input;
    std::string line;

    input.open(filename);
    if (input.fail()){
        throw(std::invalid_argument("Could not find " + filename + "."));
    }

    size_t row_idx = 0;
    size_t column_idx = 0;
    size_t idx = 0;

    data.resize(n_rows * n_columns);

    while (getline(input, line)) {
        /* Ignore empty lines or lines beginning with # */
        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        if (row_idx >= n_rows) {
            throw(std::runtime_error("Number of rows exceeds n_rows."));
        }

        double value;
        std::stringstream ss(line);

        column_idx = 0;
        idx = row_idx;

        while (ss >> value) {
            if (column_idx >= n_columns) {
                throw(
                    std::runtime_error("Number of columns exceeds n_columns."));
            }

            data.at(idx) = value;
            column_idx++;
            idx = row_idx + column_idx * n_rows;
        }
        if (column_idx != n_columns) {
            throw(std::runtime_error(
                "Number of columns does not equal n_columns."));
        }
        row_idx++;
    }
    if (row_idx != n_rows) {
        throw(std::runtime_error("Number of columns does not equal n_rows."));
    }
    input.close();
}



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
        ) 
{
    std::ofstream out(output_file);

    if (out.fail()){
        throw(std::runtime_error("Could not open " + output_file +
                                 " for writing."));
    }

    out << "# Matter power spectrum P(k) at " << n_loops << "-loop for k = " <<
        k << " (h/Mpc)\n";
    out << "#\n# Description: " << description << "\n";
    out <<    "# Git hash:    " << build_git_sha << "\n";
    out <<    "# Build time:  " << build_git_time << "\n";

    out << "#\n# Correlations computed (zero-indexed components):\n# ";
    for (auto& el : results.get_correlations()) {
        out << el <<  " ; ";
    }

    out << "\n#\n# Parameters used:\n";
    out << "# Input power spectrum read from = " << input_ps_file <<
        "\n#\n";

    out << std::scientific;
    out << "# Integration limits             = [" << q_min << "," << q_max << "]\n";
    out << "# Monte Carlo abstol, reltol     = " << cuba_epsabs << ", " <<
        cuba_epsrel << "\n";
    out << "# Monte Carlo max num. of evals  = " << cuba_maxevals << "\n";

    /* A column consists of 12 characters; 4 whitespaces in between each column */
    
    out << "#\n#" << setw(15) << "k (h/Mpc)";
    for (auto& el : results.get_correlations()) {
        out << setw(11) << "P_lin " << el;
        out << setw(5) << "P_" << n_loops << "loop " << el;
        out << setw(5) << "err_" << n_loops << "loop " << el;
    }
    out << "\n";

    out << std::setw(16) << k;
    for (size_t i = 0; i < results.get_correlations().size(); ++i) {
        out << setw(16) << results.lin_ps.at(i);
        out << setw(16) << results.non_lin_ps.at(i);
        out << setw(16) << results.errors.at(i);
    }
    out << std::endl;

    out.close();
}
