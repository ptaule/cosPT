/*
   combinatorics.hpp

   Created by Petter Taule on 08.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <ostream>
#include <vector>

class Combinations {
    private:
        int n_;
        int k_;

#if DEBUG >= 1
        /* Check that we get (n choose k) when we have gone through all
         * combinations */
        int counter = 1;
#endif

        std::vector<int> combination;
        std::vector<int> complement;

        bool next_combination();
        bool next_complement();
    public:
        Combinations(int n, int k);

        int n() const {return n_;}
        int k() const {return k_;}

        bool next();
        void reset();
        std::vector<int> get_current_combination() const {
            return combination;
        };
        std::vector<int> get_current_complement() const {
            return complement;
        }

        void write_current_combination(std::vector<int>& vec) const {
            vec = combination;
        }
        void write_current_complement(std::vector<int>& vec) const {
            vec = complement;
        }

        void rearrange_from_current(
                std::vector<int>::iterator first,
                std::vector<int>::iterator last
                ) const;
        void rearrange_from_current_combination(
                const int original[],
                int rearranged[],
                std::size_t size
                ) const;
        void rearrange_from_current_complement(
                const int original[],
                int rearranged[],
                std::size_t size
                ) const;
};


class Orderings {
    private:
        int n_;

        std::vector<Combinations> combinations_vec;
        std::vector<int> normal_ordering;
    public:
        Orderings() : n_(), combinations_vec(), normal_ordering() {}
        Orderings(int n, const std::vector<int>& group_sizes);

        Orderings& operator=(const Orderings& other) = default;
        Orderings& operator=(Orderings&& other) = default;

        int n() const {return n_;}

        bool next();
        void reset();

        void write_current(std::vector<int>& ordering) const;
        std::vector<int> get_current() const;
};


/* Variant of combinatorics problem "stars and bars", where the stars are
* unique elements (numbers). Have n elements and want to group them in m bins
* (or equivalently divide using (m-1) bars). Go through all prossibilities
* using next(), permutations within a group is not counted. */
class NumbersAndBars {
    private:
        int n_;
        int m_;

        Combinations comb;
        Orderings orderings;
        std::vector<int> current_comb;
        std::vector<int> ordering;
        std::vector<int> group_sizes;

        void next_combination_compute();
        void next_ordering_compute() { orderings.write_current(ordering); }
    public:
        NumbersAndBars(int n, int m);

        bool next();
};

std::ostream& operator<<(std::ostream& out, const Combinations& combinations);
std::ostream& operator<<(std::ostream& out, const Orderings& orderings);

#endif /* ifndef COMBINATORICS_HPP */
