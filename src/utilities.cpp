#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <sstream>

#include "../include/utilities.hpp"

using std::size_t;
using std::pow;

std::ostream& operator<<(std::ostream& out, const Pair<int>& pair) {
    out << "<" << pair.first() << "," << pair.second() << ">";
    return out;
}



std::ostream& operator<<(std::ostream& out, const Triple<int>& pair) {
    out << "<" << pair.first() << "," << pair.second() << "," << pair.third()
        << ">";
    return out;
}



#if DEBUG == 0
#define at(x) operator[](x)
#endif
/* Create config vector from label */
void label2config(int label, Vec1D<int>& config)
{
    /* Coefficients are element in {-1,0,1} */
    for (size_t i = 0; i < config.size(); ++i) {
        config.at(i) = (label % 3) - 1;
        label /= 3;
    }
}



/* Get label from config vector */
int config2label(const Vec1D<int>& config)
{
    int label = 0;
    for (size_t i = 0; i < config.size(); ++i) {
        /* Add 1 to obtain range {0,1,2} */
        label += (config.at(i) + 1) * static_cast<int>(pow(3,i));
    }
    return label;
}
#undef at



std::string label_to_string(
    int label,
    size_t n_coeffs,
    Spectrum spectrum
) {
    Vec1D<int> config(n_coeffs, 0);
    label2config(label, config);

    std::ostringstream oss;

    auto print_sign = [&](int val, const std::string& wavevector) {
        if (val == -1) oss << "-" << wavevector;
        else if (val == 1) oss << "+" << wavevector;
    };

    if (spectrum == POWERSPECTRUM) {
        print_sign(config[n_coeffs - 1], "k");
        oss << "";
        print_sign(config[n_coeffs - 2], "Q" + std::to_string(n_coeffs - 1));
    } else if (spectrum == BISPECTRUM) {
        print_sign(config[n_coeffs - 1], "ka");
        print_sign(config[n_coeffs - 2], "kb");
    }

    for (size_t i = 0; i < n_coeffs - 2; ++i) {
        print_sign(config[i], "Q" + std::to_string(i + 1));
    }
    return oss.str();
}



std::string labels_to_string(
    const int labels[],
    size_t size,
    size_t n_coeffs,
    Spectrum spectrum
    )
{
    std::ostringstream oss;
    oss << "(";
    bool first = true;
    for (size_t i = 0; i < size; ++i) {
        if (labels[i] == get_zero_label(n_coeffs)) continue;
        if (!first) oss << ", ";
        oss << label_to_string(labels[i], n_coeffs, spectrum);
        first = false;
    }
    oss << ")";
    return oss.str();
}



int get_zero_label(size_t n_coeffs) {
    const Vec1D<int> config(n_coeffs, 0);
    return config2label(config);
}



bool single_loop_label(int label, size_t n_coeffs, Spectrum spectrum)
{
    Vec1D<int> config(n_coeffs, 0);
    label2config(label, config);

    size_t current = 0;
    if (spectrum == POWERSPECTRUM) {
        /* Last label is k, which is not present when the label is pure loop
         * momentum*/
        if (config.at(n_coeffs-1) != 0) {
            return false;
        }
        current = n_coeffs - 2;
    }
    else if (spectrum == BISPECTRUM) {
        /* Two last labels are k_a and k_b, which are not present when the
         * label is pure loop momentum*/
        if (config.at(n_coeffs - 1) != 0 || config.at(n_coeffs - 2) != 0) {
            return false;
        }
        current = n_coeffs - 3;
    }

    int num_vecs_present = 0;

    /* Go through rest of coefficients and count present momenta */
    for (int i = static_cast<int>(current); i >= 0; --i) {
        if (config.at(static_cast<size_t>(i)) != 0) num_vecs_present++;
    }
    return (num_vecs_present == 1);
}



int flip_signs(int label, size_t n_coeffs)
{
    std::vector<int> config(n_coeffs);
    label2config(label, config);
    std::transform(config.begin(), config.end(), config.begin(),
                   [](int x) { return -x; });
    return config2label(config);
}
