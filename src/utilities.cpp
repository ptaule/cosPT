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
        short int label,    // in, label element in [0, N_CONFIG - 1]
        short int config[], // out, array of base 3 digits
        size_t size   // in, size of coefficients
        )
{
    // Q-coefficients are element in {-1,0,1}
    for (std::size_t i = 0; i < size - 1; ++i) {
        config[i] = (label % 3) - 1;
        label /= 3;
    }

    // k-coefficient (last coefficient) is element in {0,1}

    config[size - 1] = label % 3;
}



short int config2label(
        const short int config[], // in, array of base 3 digits
        size_t size         // in, size of coefficients
        )
{
    short int label = 0;

    // For Q-coefficients, add 1 to obtain range {0,1,2}
    for (size_t i = 0; i < size - 1; ++i) {
        label += (config[i] + 1) * pow(3,i);
    }

    // For k-coefficient, do _not_ add 1, range is {0,1}
    label += config[size-1] * pow(3,size - 1);

    return label;
}



void print_label(
        short int label, 
        short int n_coeffs, 
        std::ostream& out
        ) 
{
    short int config[N_COEFFS_MAX];
    label2config(label, config, n_coeffs);

    if (config[n_coeffs - 1] == 1) {
        out << "k";
    }
    for (int i = 0; i < n_coeffs - 1; ++i) {
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
        const short int labels[], 
        size_t size,
        short int n_coeffs, 
        short int zero_label, 
        std::ostream& out
        )
{
    out << "(";
    for (size_t i = 0; i < size; ++i) {
        if (labels[i] == zero_label) continue;
        print_label(labels[i], n_coeffs, out);
        out << ", ";
    }
    out << ")";
}



short int get_zero_label(short int n_coeffs) {
    const short int coeffs[N_COEFFS_MAX] = {0};
    return config2label(coeffs, n_coeffs);
}



bool is_fundamental(short int label, short int n_coeffs, short int n_configs) {
    // Vector is not fundamental if k is present; this is the case if
    // label >= N_CONFIGS/2
    if (label >= n_configs/2) return false;

    short int coeffs[N_COEFFS_MAX] = {0};
    label2config(label, coeffs, n_coeffs);

    short int num_vecs_present = 0;

    // The last coefficient is for k, hence we can skip this (j < N_COEFFS - 1)
    for (int i = 0; i < n_coeffs - 1; ++i) {
        if (coeffs[i] != 0) num_vecs_present++;
    }
    return (num_vecs_present == 1);
}



// Are there duplicate elements of the array? Yes, return true;
// no, return false. Do not consider elements equal to skip.
bool unique_elements(const short int array[], size_t size, short int skip) {
    for (size_t i = 0; i < size; ++i) {
        short int val = array[i];
        if (val == skip) continue;
        for (size_t j = i + 1; j < size; ++j) {
            if (array[j] == val) return false;
        }
    }
    return true;
}
