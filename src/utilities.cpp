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

void label2config(
        int label,    /* in, label element in [0, n_config - 1] */
        int config[], /* out, array of base 3 digits            */
        size_t size   /* in, size of coefficients               */
        )
{
    // Coefficients are element in {-1,0,1}
    for (size_t i = 0; i < size; ++i) {
        config[i] = (label % 3) - 1;
        label /= 3;
    }
}



int config2label(
        const int config[], // in, array of base 3 digits
        size_t size         // in, size of config[]
        )
{
    int label = 0;

    // Add 1 to obtain range {0,1,2}
    for (size_t i = 0; i < size; ++i) {
        label += (config[i] + 1) * pow(3,i);
    }
    return label;
}



void print_label(
        int label,
        int n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        )
{
    int config[N_COEFFS_MAX];
    label2config(label, config, n_coeffs);

    if (spectrum == POWERSPECTRUM) {
        if (config[n_coeffs - 1] == -1) {
            out << "-k";
        }
        else if (config[n_coeffs - 1] == 1) {
            out << "+k";
        }

        if (config[n_coeffs - 2] == -1) {
            out << "-Q" << n_coeffs - 1;
        }
        else if (config[n_coeffs - 2] == 1) {
            out << "+Q" << n_coeffs - 1;
        }
    }
    else if (spectrum == BISPECTRUM) {
        if (config[n_coeffs - 1] == -1) {
            out << "-ka";
        }
        else if (config[n_coeffs - 1] == 1) {
            out << "+ka";
        }

        if (config[n_coeffs - 2] == -1) {
            out << "-kb";
        }
        else if (config[n_coeffs - 2] == 1) {
            out << "+kb";
        }
    }

    for (int i = 0; i < n_coeffs - 2; ++i) {
        if (config[i] == 0) continue;
        else if (config[i] == -1) {
            out << "-Q" << i + 1;
        }
        else if (config[i] == 1) {
            out << "+Q" << i + 1;
        }
    }
}



void print_labels(
        const int labels[],
        size_t size,
        int n_coeffs,
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



void change_sign(int config[], std::size_t size) {
    for (size_t i = 0; i < size; ++i) {
        config[i] *= -1;
    }
}



int get_zero_label(int n_coeffs) {
    const int coeffs[N_COEFFS_MAX] = {0};
    return config2label(coeffs, n_coeffs);
}



bool single_loop_label(int label, int n_coeffs, Spectrum spectrum)
{
    int coeffs[N_COEFFS_MAX] = {0};
    label2config(label, coeffs, n_coeffs);

    int current = 0;
    if (spectrum == POWERSPECTRUM) {
        /* Last label is k, which is not present when the label is pure loop
         * momentum*/
        if (coeffs[n_coeffs-1] != 0) {
            return false;
        }
        current = n_coeffs - 2;
    }
    else if (spectrum == BISPECTRUM) {
        /* Two last labels are k_a and k_b, which are not present when the
         * label is pure loop momentum*/
        if (coeffs[n_coeffs-1] != 0 || coeffs[n_coeffs - 2] != 0) {
            return false;
        }
        current = n_coeffs - 3;
    }

    int num_vecs_present = 0;

    /* Go through rest of coefficients and count present momenta */
    for (int i = current; i >= 0; --i) {
        if (coeffs[i] != 0) num_vecs_present++;
    }
    return (num_vecs_present == 1);
}
