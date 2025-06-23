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

// Precompute binomial coefficients up to n=6
constexpr int MAX_N = 6;

constexpr std::array<std::array<int, MAX_N + 1>, MAX_N + 1> make_binomial_table() {
    std::array<std::array<int, MAX_N + 1>, MAX_N + 1> table{};
    for (std::size_t n = 0; n <= MAX_N; ++n) {
        table[n][0] = table[n][n] = 1;
        for (std::size_t k = 1; k < n; ++k) {
            table[n][k] = table[n - 1][k - 1] + table[n - 1][k];
        }
    }
    return table;
}
constexpr auto binomial_coeffs = make_binomial_table();

enum Spectrum {POWERSPECTRUM, BISPECTRUM};
enum Dynamics {EDS_SPT, EVOLVE_EDS_ICS, EVOLVE_ASYMPTOTIC_ICS};

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

// Usage: std::cout << Colors::RED << "Error" << Colors::RESET << "\n";


/* Unused parameters (to silence gcc warnings) */
#define UNUSED(x) (void)(x)

template <class T>
using Vec1D = std::vector<T>;

template <class T>
using Vec2D = std::vector<std::vector<T>>;

template <class T>
using Vec3D = std::vector<std::vector<std::vector<T>>>;

/* Strided 2D vector class. Internal flat layout for fast memory access */
template<typename T>
class Strided2DVec {
    private:
        std::size_t stride_ = 0;
        std::vector<T> data_;
    public:
        Strided2DVec() = default;
        Strided2DVec(const Strided2DVec&) = default;
        Strided2DVec(Strided2DVec&&) = default;
        Strided2DVec& operator=(const Strided2DVec&) = default;
        Strided2DVec& operator=(Strided2DVec&&) = default;
        Strided2DVec(std::size_t rows, std::size_t cols)
            : stride_(cols), data_(rows * cols, T()) {}

        void resize(std::size_t rows, std::size_t cols, const T& value = T()) {
            stride_ = cols;
            data_.resize(rows * cols, value);
        }

        std::vector<T>& data() { return data_; }
        const std::vector<T>& data() const { return data_; }
        T* pointer() { return data_.data(); }

        std::size_t stride() const { return stride_; }
        std::size_t rows() const { return data_.size() / stride_; }
        std::size_t size() const { return data_.size(); }

        // Non-const access
        auto begin() { return data_.begin(); }
        auto end()   { return data_.end(); }
        // Const access
        auto begin() const { return data_.begin(); }
        auto end()   const { return data_.end(); }
        auto cbegin() const { return data_.cbegin(); }
        auto cend()   const { return data_.cend(); }
#if DEBUG == 0
#define at(x) operator[](x)
#endif
        T& operator()(std::size_t i, std::size_t j) {
            return data_.at(i * stride_ + j);
        }
        const T& operator()(std::size_t i, std::size_t j) const {
            return data_.at(i * stride_ + j);
        }
#undef at
};


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

void label2config(int label, std::vector<int>& config);
int config2label(const std::vector<int>& config);

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
