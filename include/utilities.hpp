#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

/* Maximum number of kernel arguments at 2-loop */
#define N_KERNEL_ARGS_MAX 6

/* Number of EdS SPT components */
#define EDS_SPT_COMPONENTS 2

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

namespace Colors {
    constexpr char RESET[]   = "\x1b[0m";
    constexpr char RED[]     = "\x1b[31m";
    constexpr char GREEN[]   = "\x1b[32m";
    constexpr char YELLOW[]  = "\x1b[33m";
    constexpr char BLUE[]    = "\x1b[34m";
    constexpr char MAGENTA[] = "\x1b[35m";
    constexpr char CYAN[]    = "\x1b[36m";
}

/* Unused parameters (to silence gcc warnings) */
#define UNUSED(x) (void)(x)

/* Vector aliases */
template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;
template <class T>
using Vec3D = std::vector<std::vector<std::vector<T>>>;

template <std::size_t N>
std::array<double, N> operator+=(
    std::array<double, N>& a,
    const std::array<double, N>& b
)
{
    for (std::size_t i = 0; i < N; ++i) {
        a[i] += b[i];
    }
    return a;
}

enum Spectrum {POWERSPECTRUM, BISPECTRUM};
enum Dynamics {EDS_SPT, EVOLVE_EDS_ICS, EVOLVE_ASYMPTOTIC_ICS};

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

std::string label2string(
    int label,
    size_t n_coeffs,
    Spectrum spectrum
);

std::string labels2string(
    const int labels[],
    size_t size,
    size_t n_coeffs,
    Spectrum spectrum
);

inline std::string labels2string(
    const std::vector<int>& labels,
    size_t n_coeffs,
    Spectrum spectrum
)
{
    return labels2string(labels.data(), labels.size(), n_coeffs,
                         spectrum);
}


int get_zero_label(size_t n_coeffs);
bool single_loop_label(int label, size_t n_coeffs, Spectrum spectrum);

/* Return the label corresponding to the sign flipped configuration */
int flip_signs(int label, size_t n_coeffs);

template <typename T>
/* Are there duplicate elements of the array? Yes, return true; */
/* no, return false. Do not consider elements equal to skip. */
bool has_duplicates_excluding(const T array[], size_t size, T skip) {
    for (size_t i = 0; i < size; ++i) {
        if (array[i] == skip) continue;
        for (size_t j = i + 1; j < size; ++j) {
            if (array[j] == array[i]) return true;
        }
    }
    return false;
}

#endif /* ifndef UTILITIES_HPP */
