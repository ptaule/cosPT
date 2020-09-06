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

/* Maximum coefficient c-array size (at 2-loop) */
#define N_COEFFS_MAX 3
/* Maximum number of kernel arguments at 2-loop */
#define N_KERNEL_ARGS_MAX 5
#define SPT_COMPONENTS 2

/* Constants: */
#define PI    3.14159265359
#define TWOPI 6.28318530718

/* Macros for performant square and cubing */
#define SQUARE(x) x*x
#define CUBE(x) x*x*x

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

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

/* Vector aliases */
template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;
template <class T>
using Vec3D = std::vector<std::vector<std::vector<T>>>;

enum Spectrum {POWERSPECTRUM, BISPECTRUM};
enum Dynamics {SPT, EVOLVE_SPT_IC, EVOLVE_ASYMP_IC};

/* Utility functions */

void label2config(short int label, short int config[], std::size_t size);
short int config2label(const short int config[], std::size_t size);


void print_label(
        short int label, 
        short int n_coeffs, 
        short int zero_label, 
        std::ostream& out
        );

void print_labels(
        const short int labels[], 
        std::size_t size,
        short int n_coeffs, 
        short int zero_label, 
        std::ostream& out
        );

short int get_zero_label(short int n_coeffs);
bool is_fundamental(short int label, short int n_coeffs, short int n_configs);
bool unique_elements(const short int array[], std::size_t length, short int skip);

#endif /* ifndef UTILITIES_HPP */
