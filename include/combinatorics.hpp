/*
   combinatorics.hpp

   Created by Petter Taule on 08.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <iosfwd>
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
        Orderings(int n, const std::vector<int>& group_sizes);
        int n() const {return n_;}

        bool next();
        void reset();

        std::vector<int> get_current() const;
};

std::ostream& operator<<(std::ostream& out, const Combinations& combinations);
std::ostream& operator<<(std::ostream& out, const Orderings& orderings);

#endif /* ifndef COMBINATORICS_HPP */
