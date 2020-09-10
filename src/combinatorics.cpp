/*
   combinatorics.cpp

   Created by Petter Taule on 08.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <ostream>
#include <vector>
#include <numeric>
#include <stdexcept>

#include <gsl/gsl_sf.h>

#include "../include/combinatorics.hpp"

template <class T>
using Vec1D = std::vector<T>;

using std::size_t;

Combinations::Combinations(int n, int k)
    : n(n), k(k)
{
    if (k < 0 || k > n) {
        throw(std::invalid_argument("Combinations::Combinations(): k >= 0 && k <= n is required."));
    }
    combination.resize(k);
    complement.resize(n - k);
    std::iota(combination.begin(), combination.end(), 0);
    std::iota(complement.begin(), complement.end(), k);
}



bool Combinations::next_combination() {
    if (combination.empty()) {
        return false;
    }

    /* Last element */
    if (combination.back() < (n - 1)) {
        ++combination.back();
        return true;
    }
    /* Go through the rest */
    for (int i = k - 2; i >= 0; --i) {
        /* Element i cannot be greater than (n - (k - i)) */
        if (combination[i] < (n - k + i)) {
            /* Increment this element, elements right of this is set in
             * increasing order */
            std::iota(combination.begin() + i, combination.end(), ++combination[i]);
            return true;
        }
    }
    /* We have gone through all combinations */
    return false;
}



bool Combinations::next_complement() {
    if (complement.empty()) {
        return false;
    }
    /* Go through complement backwards, if number is not 0 or one more
     * than previous number, decrement and return */

    for (int i = complement.size() - 1; i > 0; --i) {
        if (complement[i] > (complement[i-1] + 1)) {
            /* Decrement this element, elements right of this is set to their
             * maximal value */
            --complement[i];
            std::iota(complement.begin() + i + 1, complement.end(), k + i + 1);
            return true;
        }
    }
    /* First element */
    if (complement[0] > 0) {
        --complement[0];
        std::iota(complement.begin() + 1, complement.end(), k + 1);
        return true;
    }
    /* We have gone through all complements */
    return false;
}



bool Combinations::next() {
    bool next_comb = next_combination();
    bool next_comp = next_complement();

    /* If bool variables differ, throw error */
    if (next_comb ^ next_comp) {
        throw(std::logic_error(
            "Combinations::next(): inconsistency between combination and complement."));
    }
#if DEBUG >= 1
    if (next_comb) {
        ++counter;
    }
    else {
        if (counter != gsl_sf_choose(n,k)) {
            throw(std::logic_error(
                "Combinations::next(): Did not create (n choose k) combinations."));
        }
    }
#endif
    return next_comb;
}



void Combinations::reset()
{
    std::iota(combination.begin(), combination.end(), 0);
    std::iota(complement.begin(), complement.end(), k);
}



void Combinations::rearrange_from_current(
        Vec1D<short int>::iterator first,
        Vec1D<short int>::iterator last
        ) const
{
    if (first == last) {
        throw(std::logic_error(
                    "Combinations::rearrange_from_current(): first == last."));
    }
    Vec1D<short int> copy(first, last);

    if (copy.size() != static_cast<size_t>(n)) {
        throw(std::logic_error("Combinations::rearrange_from_current(): size of subvector from first,last is not n."));
    }

    for (int i = 0; i < k; ++i) {
        *(first + i) = copy[combination[i]];
    }
    for (int i = k; i < n; ++i) {
        *(first + i) = copy[complement[i-k]];
    }
}



void Combinations::rearrange_from_current_combination(
        const short int original[],
        short int rearranged[],
        size_t size
        ) const
{
    if (size != static_cast<size_t>(k)) {
        throw(std::invalid_argument("Combinations::rearrange_from_current_combination(): size != k."));
    }
    for (int i = 0; i < k; ++i) {
        rearranged[i] = original[combination[i]];
    }
}



void Combinations::rearrange_from_current_complement(
        const short int original[],
        short int rearranged[],
        size_t size
        ) const
{
    if (size != static_cast<size_t>(n - k)) {
        throw(std::invalid_argument("Combinations::rearrange_from_current_complement(): size != k."));
    }
    for (int i = 0; i < n - k; ++i) {
        rearranged[i] = original[complement[i]];
    }
}



Orderings::Orderings(int n, const Vec1D<short int>& group_sizes)
    : n(n)
{
    if (n < 1) {
        throw(std::invalid_argument("Orderings::Orderings(): n must be greater than 0."));
    }

    int sum = 0;
    for (auto& el : group_sizes) {
        if (el <= 0) {
            throw(std::invalid_argument("Orderings::Orderings(): group_size must be strictly positive."));
        }
        sum += el;
    }
    if ((sum != n)) {
        throw(std::invalid_argument("Orderings::Orderings(): sum of group sizes does not equal n."));
    }

    /* Start with increasing order (0,...,n) */
    normal_ordering.resize(n);
    std::iota(normal_ordering.begin(), normal_ordering.end(), 0);

    combinations_vec.reserve(group_sizes.size());

    int remaining = n;
    for (size_t i = 0; i < group_sizes.size(); ++i) {
        combinations_vec.push_back(Combinations(remaining, group_sizes[i]));
        remaining -= group_sizes[i];
    }
}



bool Orderings::next()
{
    /* Call next on innermost combination */
    int index = combinations_vec.size() - 1;
    while (index > 0) {
        if (combinations_vec.at(index).next()) {
            return true;
        }
        else {
            /* Reset the innermore combination and move to the outermore */
            combinations_vec.at(index).reset();
            --index;
        }
    }
    /* Outermost combination: if next is false we have gone through all
     * orderings and we return false */
    if (combinations_vec.at(0).next()) {
        return true;
    }
    else {
        return false;
    }
}



void Orderings::reset()
{
    for (auto& el : combinations_vec) {
        el.reset();
    }
}



std::vector<short int> Orderings::get_current() const
{
    Vec1D<short int> ordering(normal_ordering);

    int cursor = 0;
    for (size_t i = 0; i < combinations_vec.size(); ++i) {
        combinations_vec.at(i).rearrange_from_current(
                ordering.begin() + cursor, ordering.end());
        cursor += combinations_vec.at(i).get_k();
    }
    return ordering;
}



std::ostream& operator<<(std::ostream& out, const Combinations& combinations) {
    Vec1D<short int> combination = combinations.get_current_combination();
    Vec1D<short int> complement = combinations.get_current_complement();

    for (auto& el : combination) {
        out << el << "  ";
    }
    out << "\t | \t";
    for (auto& el : complement) {
        out << el << "  ";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const Orderings& orderings) {
    Vec1D<short int> ordering = orderings.get_current();

    for (auto& el : ordering) {
       out << el << "  ";
    }
    return out;
}
