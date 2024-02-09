#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../include/utilities.hpp"
#include "../include/parameters.hpp"
#include "../include/io.hpp"

using std::size_t;
using std::string;

/* Read delimited (with spaces/tabs) file and store values in Vec2D<double>
 * data (row-by-row, any existing content removed). Skip empty lines or lines
 * beginning with #. */
void read_delimited_file(
        const string& filename, /* file to read */
        Vec2D<double>& data     /* out: data */
        )
{
    std::ifstream input(filename);
    string line;

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
    std::ofstream out(cfg.get<string>("output_file"));

    if (out.fail()){
        throw(std::runtime_error("Could not open \"" + cfg.get<string>("output_file") +
                                 "\" for writing."));
    }

    out << cfg;

    /* A column consists of 14 characters; 4 whitespaces in between each column */
    if (cfg.get<Spectrum>("spectrum") == POWERSPECTRUM) {
        out << "#" << setw(17) << "k";

        if (cfg.get<bool>("rsd")) {
            for (int i = 0; i < 3; ++i) {
                out << setw(14) << "P_l" << 2*i << "_L0";
                if (cfg.get<int>("n_loops") > 0) {
                    out << setw(14)  << "P_l" << 2*i << "_L"   << cfg.get<int>("n_loops");
                    out << setw(14)  << "err_l" << 2*i << "_L" << cfg.get<int>("n_loops");
                }
            }
        }
        else {
            for (auto& el : cfg.pair_correlations()) {
                out << setw(13) << "P_L0 " << el;
                if (cfg.get<int>("n_loops") > 0) {
                    out << setw(11)  << "P_L"   << cfg.get<int>("n_loops") << " " << el;
                    out << setw(11)  << "err_L" << cfg.get<int>("n_loops") << " " << el;
                }
            }
        }
    }
    else {
        out << "#" << setw(17) << "k_a";
        out << setw(18) << "k_b";
        out << setw(18) << "k_c";

        for (auto& el : cfg.triple_correlations()) {
            out << setw(13) << "B_L0 " << el;
            if (cfg.get<int>("n_loops") > 0) {
                out << setw(13)  << "B_L"   << cfg.get<int>("n_loops") << " " << el;
                out << setw(13)  << "err_L" << cfg.get<int>("n_loops") << " " << el;
            }
        }
    }

    out << "\n";
    out << std::setw(18) << cfg.get<double>("k_a");

    if (cfg.get<Spectrum>("spectrum") == POWERSPECTRUM) {
        if (cfg.get<bool>("rsd")) {
            for (size_t i = 0; i < 3; ++i) {
                out << setw(18) << tree_level_result.at(i);
                if (cfg.get<int>("n_loops") > 0) {
                    out << setw(18) << loop_result.at(i);
                    out << setw(18) << errors.at(i);
                }
            }
        }
        else {
            for (size_t i = 0; i < cfg.pair_correlations().size(); ++i) {
                out << setw(18) << tree_level_result.at(i);
                if (cfg.get<int>("n_loops") > 0) {
                    out << setw(18) << loop_result.at(i);
                    out << setw(18) << errors.at(i);
                }
            }
        }
    }
    else {
        out << setw(18) << cfg.get<double>("k_b");
        out << setw(18) << cfg.get<double>("k_c");

        for (size_t i = 0; i < cfg.triple_correlations().size(); ++i) {
            out << setw(18) << tree_level_result.at(i);
            if (cfg.get<int>("n_loops") > 0) {
                out << setw(18) << loop_result.at(i);
                out << setw(18) << errors.at(i);
            }
        }
    }
    out << std::endl;

    out.close();
}
