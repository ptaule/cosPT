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
        int n;
        int k;

#if DEBUG >= 1
        /* Check that we get (n choose k) when we have gone through all
         * combinations */
        int counter = 1;
#endif

        std::vector<short int> combination;
        std::vector<short int> complement;

        bool next_combination();
        bool next_complement();
    public:
        Combinations(int n, int k);

        int get_n() const {return n;}
        int get_k() const {return k;}

        bool next();
        void reset();
        std::vector<short int> get_current_combination() const {
            return combination;
        };
        std::vector<short int> get_current_complement() const {
            return complement;
        }

        void rearrange_from_current(
                std::vector<short int>::iterator first,
                std::vector<short int>::iterator last
                ) const;
        void rearrange_from_current_combination(
                const short int original[],
                short int rearranged[],
                std::size_t size
                ) const;
        void rearrange_from_current_complement(
                const short int original[],
                short int rearranged[],
                std::size_t size
                ) const;
};

class Orderings {
    private:
        int n;

        std::vector<Combinations> combinations_vec;
        std::vector<short int> normal_ordering;
    public:
        Orderings(int n, const std::vector<short int>& group_sizes);

        bool next();
        void reset();

        std::vector<short int> get_current() const;
};

std::ostream& operator<<(std::ostream& out, const Combinations& combinations);
std::ostream& operator<<(std::ostream& out, const Orderings& orderings);

#endif /* ifndef COMBINATORICS_HPP */
