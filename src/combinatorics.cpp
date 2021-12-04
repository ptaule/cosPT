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

/* Turn off vector bounds check if not in debug-mode */
#if DEBUG == 0
#define at(x) operator[](x)
#endif

template <class T>
using Vec1D = std::vector<T>;

using std::size_t;

Combinations::Combinations(int n, int k)
    : n_(n), k_(k)
{
    if (k_ < 0 || k_ > n) {
        throw(std::invalid_argument(
            "Combinations::Combinations(): k >= 0 && k <= n is required."));
    }
    combination.resize(static_cast<size_t>(k_));
    complement.resize(static_cast<size_t>(n_ - k_));
    std::iota(combination.begin(), combination.end(), 0);
    std::iota(complement.begin(), complement.end(), k_);
}



bool Combinations::next_combination() {
    if (combination.empty()) {
        return false;
    }

    /* Last element */
    if (combination.back() < (n_ - 1)) {
        ++combination.back();
        return true;
    }
    /* Go through the rest */
    for (int i = k_ - 2; i >= 0; --i) {
        /* Element i cannot be greater than (n - (k - i)) */
        if (combination.at(static_cast<size_t>(i)) < (n_ - k_ + i)) {
            /* Increment this element, elements right of this is set in
             * increasing order */
            std::iota(combination.begin() + i, combination.end(),
                    ++combination.at(static_cast<size_t>(i)));
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

    for (size_t i = complement.size() - 1; i > 0; --i) {
        if (complement.at(i) > (complement.at(i-1) + 1)) {
            /* Decrement this element, elements right of this is set to their
             * maximal value */
            --complement.at(i);
            std::iota(complement.begin() + static_cast<int>(i) + 1,
                    complement.end(), k_ + static_cast<int>(i) + 1);
            return true;
        }
    }
    /* First element */
    if (complement.at(0) > 0) {
        --complement.at(0);
        std::iota(complement.begin() + 1, complement.end(), k_ + 1);
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
        throw(std::logic_error("Combinations::next(): inconsistency between "
                               "combination and complement."));
    }
#if DEBUG >= 1
    if (next_comb) {
        ++counter;
    }
    else {
        if (counter != gsl_sf_choose(
                    static_cast<unsigned int>(n_),
                    static_cast<unsigned int>(k_)
                    )
                ) {
            throw(std::logic_error("Combinations::next(): Did not create (n "
                                   "choose k) combinations."));
        }
    }
#endif
    return next_comb;
}



void Combinations::reset()
{
    std::iota(combination.begin(), combination.end(), 0);
    std::iota(complement.begin(), complement.end(), k_);
#if DEBUG >= 1
    counter = 1;
#endif
}



void Combinations::rearrange_from_current(
        Vec1D<int>::iterator first,
        Vec1D<int>::iterator last
        ) const
{
    if (first == last) {
        throw(std::logic_error(
            "Combinations::rearrange_from_current(): first == last."));
    }
    Vec1D<int> copy(first, last);

    if (copy.size() != static_cast<size_t>(n_)) {
        throw(std::logic_error("Combinations::rearrange_from_current(): size "
                               "of subvector from first,last is not n."));
    }

    for (int i = 0; i < k_; ++i) {
        *(first + i) =
            copy.at(static_cast<size_t>(combination.at(static_cast<size_t>(i))));
    }
    for (int i = k_; i < n_; ++i) {
        *(first + i) =
            copy.at(static_cast<size_t>(complement.at(static_cast<size_t>(i - k_))));
    }
}



void Combinations::rearrange_from_current_combination(
        const int original[],
        int rearranged[],
        size_t size
        ) const
{
    if (size != static_cast<size_t>(k_)) {
        throw(std::invalid_argument(
            "Combinations::rearrange_from_current_combination(): size != k."));
    }
    for (size_t i = 0; i < static_cast<size_t>(k_); ++i) {
        rearranged[i] = original[combination.at(i)];
    }
}



void Combinations::rearrange_from_current_complement(
        const int original[],
        int rearranged[],
        size_t size
        ) const
{
    if (size != static_cast<size_t>(n_ - k_)) {
        throw(std::invalid_argument(
            "Combinations::rearrange_from_current_complement(): size != k."));
    }
    for (size_t i = 0; i < static_cast<size_t>(n_ - k_); ++i) {
        rearranged[i] = original[complement.at(i)];
    }
}



/* Orderings code is not performance critical, therefore always check bounds */
#undef at



Orderings::Orderings(int n, const Vec1D<int>& group_sizes)
    : n_(n)
{
    if (n < 1) {
        throw(std::invalid_argument(
            "Orderings::Orderings(): n must be greater than 0."));
    }

    int sum = 0;
    for (auto& el : group_sizes) {
        if (el <= 0) {
            throw(std::invalid_argument("Orderings::Orderings(): group_size "
                                        "must be strictly positive."));
        }
        sum += el;
    }
    if ((sum != n)) {
        throw(std::invalid_argument(
            "Orderings::Orderings(): sum of group sizes does not equal n."));
    }

    /* Start with increasing order (0,...,n) */
    normal_ordering.resize(static_cast<size_t>(n));
    std::iota(normal_ordering.begin(), normal_ordering.end(), 0);

    combinations_vec.reserve(group_sizes.size());

    int remaining = n;
    for (size_t i = 0; i < group_sizes.size(); ++i) {
        combinations_vec.push_back(Combinations(remaining, group_sizes.at(i)));
        remaining -= group_sizes.at(i);
    }
}



bool Orderings::next()
{
    /* Call next on innermost combination */
    size_t index = combinations_vec.size() - 1;
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



void Orderings::write_current(Vec1D<int>& ordering) const {
    ordering = normal_ordering;

    int cursor = 0;
    for (size_t i = 0; i < combinations_vec.size(); ++i) {
        combinations_vec.at(i).rearrange_from_current(
                ordering.begin() + cursor, ordering.end());
        cursor += combinations_vec.at(i).k();
    }
}



Vec1D<int> Orderings::get_current() const
{
    Vec1D<int> ordering(n_);
    write_current(ordering);
    return ordering;
}



std::ostream& operator<<(std::ostream& out, const Combinations& combinations) {
    Vec1D<int> combination = combinations.get_current_combination();
    Vec1D<int> complement = combinations.get_current_complement();

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
    Vec1D<int> ordering = orderings.get_current();

    for (auto& el : ordering) {
       out << el << "  ";
    }
    return out;
}
