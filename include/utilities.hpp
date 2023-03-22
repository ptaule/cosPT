/*
   utilities.hpp

   Created by Petter Taule on 28.08.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cstddef>
#include <vector>
#include <iosfwd>

/* Maximum number of kernel arguments at 2-loop */
#define N_KERNEL_ARGS_MAX 6

/* Number of components */
#define EDS_SPT_COMPONENTS 2
#define COMPONENTS 4

/* Constants: */
#define SQRT2 1.41421356237
#define PI    3.14159265359
#define TWOPI 6.28318530718

/* Macros for performant powers */
#define SQUARE(x) (x)*(x)
#define CUBE(x) (x)*(x)*(x)
#define POW4(x) (x)*(x)*(x)*(x)

/* Debug modes (if not set by compile options): */
/* 1: Perform additional checks during runtime. */
/* 2: In addition, print info during */
/* runtime.  */
#ifndef DEBUG
#define DEBUG 0
#endif

// Which GSL ODE routine to use
#define ODE_ROUTINE gsl_odeiv2_step_rkf45

// Various colors
#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"
#define COLOR_RESET   "\x1b[0m"

/* Unused parameters (to silence gcc warnings) */
#define UNUSED(x) (void)(x)

/* Vector aliases */
template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;
template <class T>
using Vec3D = std::vector<std::vector<std::vector<T>>>;

enum Spectrum {POWERSPECTRUM, BISPECTRUM};
enum Dynamics {EDS_SPT, EVOLVE_IC_EDS, EVOLVE_IC_ASYMP};

/* Pair and triple classes */
template <typename T>
class Pair {
    private:
        T a_;
        T b_;
    public:
        Pair() = default;
        Pair(const Pair&) = default;
        Pair(T a, T b) : a_(a), b_(b) {}

        const T& a() const {return a_;}
        const T& b() const {return b_;}

        T& a() {return a_;}
        T& b() {return b_;}

        const T& first() const {return a_;}
        const T& second() const {return b_;}

        T& first() {return a_;}
        T& second() {return b_;}
};


template <typename T>
class Triple {
    private:
        T a_;
        T b_;
        T c_;
    public:
        Triple() = default;
        Triple(const Triple&) = default;
        Triple(T a, T b, T c) : a_(a), b_(b), c_(c) {}

        const T& a() const {return a_;}
        const T& b() const {return b_;}
        const T& c() const {return c_;}

        T& a() {return a_;}
        T& b() {return b_;}
        T& c() {return c_;}

        const T& first() const {return a_;}
        const T& second() const {return b_;}
        const T& third() const {return c_;}

        T& first() {return a_;}
        T& second() {return b_;}
        T& third() {return c_;}
};

std::ostream& operator<<(std::ostream& out, const Pair<int>& pair);
std::ostream& operator<<(std::ostream& out, const Triple<int>& pair);

/* Utility functions */

void label2config(int label, Vec1D<int>& config);
int config2label(const Vec1D<int>& config);

void print_label(
        int label,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        );

void print_labels(
        const int labels[],
        size_t size,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        );

void print_labels(
        const Vec1D<int>& labels,
        size_t n_coeffs,
        Spectrum spectrum,
        std::ostream& out
        );


int get_zero_label(size_t n_coeffs);
bool single_loop_label(int label, size_t n_coeffs, Spectrum spectrum);

template <typename T>
void flip_signs(const Vec1D<T>& input, Vec1D<T>& result) {
    for (std::size_t i = 0; i < input.size(); ++i) {
        result.at(i) = -input.at(i);
    }
}

/* Return the label corresponding to the sign flipped configuration */
int flip_signs(int label, size_t n_coeffs);

template <typename T>
/* Are there duplicate elements of the array? Yes, return true; */
/* no, return false. Do not consider elements equal to skip. */
bool unique_elements(const T array[], size_t size, T skip) {
    for (size_t i = 0; i < size; ++i) {
        T val = array[i];
        if (val == skip) continue;
        for (size_t j = i + 1; j < size; ++j) {
            if (array[j] == val) return false;
        }
    }
    return true;
}

#endif /* ifndef UTILITIES_HPP */
