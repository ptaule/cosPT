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


using std::setw;

void write_results(
        const Config& cfg,
        const Vec1D<double>& lin_ps,
        const Vec1D<double>& non_lin_ps,
        const Vec1D<double>& errors
        )
{
    std::ofstream out(cfg.output_file());

    if (out.fail()){
        throw(std::runtime_error("Could not open " + cfg.output_file() +
                                 " for writing."));
    }

    out << cfg;

    /* A column consists of 14 characters; 4 whitespaces in between each column */

    if (cfg.spectrum() == POWERSPECTRUM) {
        out << "#\n#" << setw(17) << "k (h/Mpc)";
        for (auto& el : cfg.pair_correlations()) {
            out << setw(13) << "P_lin " << el;
            out << setw(7)  << "P_"     << cfg.n_loops() << "loop " << el;
            out << setw(7)  << "err_"   << cfg.n_loops() << "loop " << el;
        }
    }
    else {
        out << "#\n#" << setw(17) << "k_a (h/Mpc)";
        out << setw(18) << "k_b (h/Mpc)";
        out << setw(18) << "cos_ab";

        for (auto& el : cfg.triple_correlations()) {
            out << setw(13) << "B_lin " << el;
            out << setw(7)  << "B_"     << cfg.n_loops() << "loop " << el;
            out << setw(7)  << "err_"   << cfg.n_loops() << "loop " << el;
        }
    }

    out << "\n";
    out << std::setw(18) << cfg.k_a();

    if (cfg.spectrum() == POWERSPECTRUM) {
        for (size_t i = 0; i < cfg.pair_correlations().size(); ++i) {
            out << setw(18) << lin_ps.at(i);
            out << setw(18) << non_lin_ps.at(i);
            out << setw(18) << errors.at(i);
        }
    }
    else {
        out << setw(18) << cfg.k_b();
        out << setw(18) << cfg.cos_ab();

        for (size_t i = 0; i < cfg.triple_correlations().size(); ++i) {
            out << setw(18) << lin_ps.at(i);
            out << setw(18) << non_lin_ps.at(i);
            out << setw(18) << errors.at(i);
        }
    }
    out << std::endl;

    out.close();
}
