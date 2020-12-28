/*
   utilities.cpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <cstddef>
#include <ostream>
#include <cmath>

#include "../include/utilities.hpp"

using std::size_t;
using std::pow;

std::ostream& operator<<(std::ostream& out, const Pair<int>& pair) {
    out << "<" << pair.first() << "," << pair.second() << ">";
    return out;
}



std::ostream& operator<<(std::ostream& out, const Triple<int>& pair) {
    out << "<" << pair.first() << "," << pair.second() << "," << pair.third()
        << ">";
    return out;
}



#if DEBUG == 0
#define at(x) operator[](x)
#endif
/* Create config vector from label */
void label2config(int label, Vec1D<int>& config)
{
    /* Coefficients are element in {-1,0,1} */
    for (size_t i = 0; i < config.size(); ++i) {
        config.at(i) = (label % 3) - 1;
        label /= 3;
    }
}



/* Get label from config vector */
int config2label(const Vec1D<int>& config)
{
    int label = 0;
    for (size_t i = 0; i < config.size(); ++i) {
        /* Add 1 to obtain range {0,1,2} */
        label += (config.at(i) + 1) * static_cast<int>(pow(3,i));
    }
    return label;
}
#undef at



void print_label(
        int label,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        )
{
    Vec1D<int> config(n_coeffs, 0);
    label2config(label, config);

    if (spectrum == POWERSPECTRUM) {
        if (config.at(n_coeffs - 1) == -1) {
            out << "-k";
        }
        else if (config.at(n_coeffs - 1) == 1) {
            out << "+k";
        }

        if (config.at(n_coeffs - 2) == -1) {
            out << "-Q" << n_coeffs - 1;
        }
        else if (config.at(n_coeffs - 2) == 1) {
            out << "+Q" << n_coeffs - 1;
        }
    }
    else if (spectrum == BISPECTRUM) {
        if (config.at(n_coeffs - 1) == -1) {
            out << "-ka";
        }
        else if (config.at(n_coeffs - 1) == 1) {
            out << "+ka";
        }

        if (config.at(n_coeffs - 2) == -1) {
            out << "-kb";
        }
        else if (config.at(n_coeffs - 2) == 1) {
            out << "+kb";
        }
    }

    for (size_t i = 0; i < n_coeffs - 2; ++i) {
        if (config.at(i) == 0) continue;
        else if (config.at(i) == -1) {
            out << "-Q" << i + 1;
        }
        else if (config.at(i) == 1) {
            out << "+Q" << i + 1;
        }
    }
}



void print_labels(
        const int labels[],
        size_t size,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        )
{
    out << "(";
    for (size_t i = 0; i < size; ++i) {
        if (labels[i] == get_zero_label(n_coeffs)) {
            continue;
        }
        print_label(labels[i], n_coeffs, spectrum, out);
        out << ", ";
    }
    out << ")";
}



void print_labels(
        const Vec1D<int>& labels,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        )
{
    return print_labels(labels.data(), labels.size(), n_coeffs, spectrum, out);
}



int get_zero_label(size_t n_coeffs) {
    const Vec1D<int> config(n_coeffs, 0);
    return config2label(config);
}



bool single_loop_label(int label, size_t n_coeffs, Spectrum spectrum)
{
    Vec1D<int> config(n_coeffs, 0);
    label2config(label, config);

    size_t current = 0;
    if (spectrum == POWERSPECTRUM) {
        /* Last label is k, which is not present when the label is pure loop
         * momentum*/
        if (config.at(n_coeffs-1) != 0) {
            return false;
        }
        current = n_coeffs - 2;
    }
    else if (spectrum == BISPECTRUM) {
        /* Two last labels are k_a and k_b, which are not present when the
         * label is pure loop momentum*/
        if (config.at(n_coeffs - 1) != 0 || config.at(n_coeffs - 2) != 0) {
            return false;
        }
        current = n_coeffs - 3;
    }

    int num_vecs_present = 0;

    /* Go through rest of coefficients and count present momenta */
    for (int i = static_cast<int>(current); i >= 0; --i) {
        if (config.at(static_cast<size_t>(i)) != 0) num_vecs_present++;
    }
    return (num_vecs_present == 1);
}



int flip_signs(int label, size_t n_coeffs) {
    Vec1D<int> config(n_coeffs);
    label2config(label, config);
    flip_signs(config, config);
    return config2label(config);
}
