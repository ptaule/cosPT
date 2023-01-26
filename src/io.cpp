/*
   io.cpp

   Created by Petter Taule on 03.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../include/utilities.hpp"
#include "../include/parameters.hpp"
#include "../include/version.hpp"
#include "../include/io.hpp"

using std::size_t;

/* Read delimited (with spaces/tabs) file and store values in Vec2D<double>
 * data (row-by-row, any existing content removed). Skip empty lines or lines
 * beginning with #. */
void read_delimited_file(
        const std::string& filename, /* file to read */
        Vec2D<double>& data          /* out: data */
        )
{
    std::ifstream input(filename);
    std::string line;

    if (input.fail()){
        throw(std::invalid_argument("Could not find \"" + filename + "\"."));
    }

    /* Clear any existing data */
    data.clear();

    /* Row counter */
    size_t i = 0;

    /* Step through rows */
    while (getline(input, line)) {
        /* Ignore empty lines or lines beginning with # */
        if (line.empty() || line.at(0) == '#') {
            continue;
        }

        data.push_back(Vec1D<double>());

        std::stringstream ss(line);
        double value;

        while (ss >> value) {
            data.at(i).push_back(value);
        }
        ++i;
    }
    input.close();
}



using std::setw;

void write_results(
        const Config& cfg,
        const Vec1D<double>& tree_level_result,
        const Vec1D<double>& loop_result,
        const Vec1D<double>& errors
        )
{
    std::ofstream out(cfg.output_file());

    if (out.fail()){
        throw(std::runtime_error("Could not open \"" + cfg.output_file() +
                                 "\" for writing."));
    }

    out << cfg;

    /* A column consists of 14 characters; 4 whitespaces in between each column */

    if (cfg.spectrum() == POWERSPECTRUM) {
        out << "#\n#" << setw(17) << "k";

        if (cfg.rsd()) {
            for (int i = 0; i < 3; ++i) {
                out << setw(14) << "P_l" << 2*i << "_L0";
                if (cfg.n_loops() > 0) {
                    out << setw(14)  << "P_l" << 2*i << "_L"   << cfg.n_loops();
                    out << setw(14)  << "err_l" << 2*i << "_L" << cfg.n_loops();
                }
            }
        }
        else {
            for (auto& el : cfg.pair_correlations()) {
                out << setw(13) << "P_L0 " << el;
                if (cfg.n_loops() > 0) {
                    out << setw(13)  << "P_L"   << cfg.n_loops() << " " << el;
                    out << setw(13)  << "err_L" << cfg.n_loops() << " " << el;
                }
            }
        }
    }
    else {
        out << "#\n#" << setw(17) << "k_a";
        out << setw(18) << "k_b";
        out << setw(18) << "k_c";

        for (auto& el : cfg.triple_correlations()) {
            out << setw(13) << "B_L0 " << el;
            if (cfg.n_loops() > 0) {
                out << setw(13)  << "B_L"   << cfg.n_loops() << " " << el;
                out << setw(13)  << "err_L" << cfg.n_loops() << " " << el;
            }
        }
    }

    out << "\n";
    out << std::setw(18) << cfg.k_a();

    if (cfg.spectrum() == POWERSPECTRUM) {
        if (cfg.rsd()) {
            for (size_t i = 0; i < 3; ++i) {
                out << setw(18) << tree_level_result.at(i);
                if (cfg.n_loops() > 0) {
                    out << setw(18) << loop_result.at(i);
                    out << setw(18) << errors.at(i);
                }
            }
        }
        else {
            for (size_t i = 0; i < cfg.pair_correlations().size(); ++i) {
                out << setw(18) << tree_level_result.at(i);
                if (cfg.n_loops() > 0) {
                    out << setw(18) << loop_result.at(i);
                    out << setw(18) << errors.at(i);
                }
            }
        }
    }
    else {
        out << setw(18) << cfg.k_b();
        out << setw(18) << cfg.k_c();

        for (size_t i = 0; i < cfg.triple_correlations().size(); ++i) {
            out << setw(18) << tree_level_result.at(i);
            if (cfg.n_loops() > 0) {
                out << setw(18) << loop_result.at(i);
                out << setw(18) << errors.at(i);
            }
        }
    }
    out << std::endl;

    out.close();
}
