#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <stdexcept>
#include <string>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

/* Wrapper classes around GSL 1d and 2d spline
 * interpolation */

class Interpolation1D {
    private:
        double x_min = 0;
        double x_max = 0;

        gsl_spline* spline;
        gsl_interp_accel* acc;
        const gsl_interp_type* type;

        void ensure_increasing_x(
                std::vector<double>& x,
                std::vector<double>& y
                );

        void initialize(
                std::vector<double> x,
                std::vector<double> y,
                double factor
                );
    public:
        Interpolation1D() : spline(nullptr), acc(nullptr), type(nullptr) {}

        /* Optional arguments:
         * - factor which multiplies y-values before interpolation
         * - gsl_interp_type: GSL interpolation type */
        Interpolation1D(
                const std::vector<double>& x,
                const std::vector<double>& y,
                double factor = 1,
                const gsl_interp_type* type = gsl_interp_cspline
                );
        Interpolation1D(
                const std::string& filename,
                double factor = 1,
                const gsl_interp_type* type = gsl_interp_cspline
                );
        Interpolation1D(const Interpolation1D&) = delete;
        Interpolation1D& operator=(const Interpolation1D&) = delete;
        Interpolation1D(Interpolation1D&& other);
        Interpolation1D& operator=(Interpolation1D&& other);

        double x_minimum() const { return x_min; }
        double x_maximum() const { return x_max; }

        /* Evaluate spline */
        double operator()(double x) const {
#if DEBUG == 0
            return gsl_spline_eval(spline, x, acc);
#else
            /* Check limits of interpolation region */
            if (x < x_min || x > x_max) {
                throw std::runtime_error(
                    "Interpolation1D::operator(): x = " + std::to_string(x) +
                    " is outside interpolation area (" + std::to_string(x_min) +
                    ", " + std::to_string(x_max) + ").");
            }
            return gsl_spline_eval(spline, x, acc);
#endif
        }

        /* Evaluate spline when min < x < max or else return 0. This is e.g.
         * useful for the power spectrum spline, calling e.g. P(k-q) can
         * evaluate the power spectrum outside the limits of the q-integration.
         * */
        double operator()(double x, double min, double max) const {
            if (x > max || x < min) {
                return 0;
            }
            else {
#if DEBUG == 0
                return gsl_spline_eval(spline, x, acc);
#else
                /* Check limits of interpolation region */
                if (x < x_min || x > x_max) {
                    throw std::runtime_error(
                        "Interpolation1D::operator(): x = " + std::to_string(x) +
                        " is outside interpolation area (" +
                        std::to_string(x_min) + ", " + std::to_string(x_max) +
                        ").");
                }
                return gsl_spline_eval(spline, x, acc);
#endif
            }
        }

        ~Interpolation1D() {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        }
};


class Interpolation2D {
    private:
        double x_min = 0;
        double x_max = 0;
        double y_min = 0;
        double y_max = 0;

        gsl_spline2d* spline;
        gsl_interp_accel* x_acc;
        gsl_interp_accel* y_acc;
        const gsl_interp2d_type* type;

        void ensure_increasing_x_y(
                std::vector<double>& x,
                std::vector<double>& y,
                std::vector<std::vector<double>>& z
                );

        void initialize(
                std::vector<double> x,
                std::vector<double> y,
                std::vector<std::vector<double>> z,
                double factor
                );

        std::size_t idx2d_to_idx1d(
                std::size_t x_idx,
                std::size_t x_size,
                std::size_t y_idx
                )
        { return x_idx + y_idx*x_size; }

    public:
        Interpolation2D() :
            spline(nullptr), x_acc(nullptr), y_acc(nullptr), type(nullptr) {}
        Interpolation2D(
                const std::vector<double>& x,
                const std::vector<double>& y,
                const std::vector<std::vector<double>>& z,
                double factor = 1,
                const gsl_interp2d_type* type = gsl_interp2d_bicubic
                );
        Interpolation2D(
                const std::string& x_grid_file,
                const std::string& y_grid_file,
                const std::string& data_file,
                double factor = 1,
                const gsl_interp2d_type* type = gsl_interp2d_bicubic
                );
        Interpolation2D(const Interpolation2D&) = delete;
        Interpolation2D& operator=(const Interpolation2D&) = delete;
        Interpolation2D(Interpolation2D&&) noexcept;
        Interpolation2D& operator=(Interpolation2D&&);

        /* Evaluate spline */
        double operator()(double x, double y) const {
#if DEBUG == 0
            return gsl_spline2d_eval(spline, x, y, x_acc, y_acc);
#else
            /* Check limits of interpolation region */
            if (x < x_min || x > x_max) {
                throw std::runtime_error(
                    "Interpolation2D::operator():x = " + std::to_string(x) +
                    " is outside interpolation area (" + std::to_string(x_min) +
                    ", " + std::to_string(x_max) + ").");
            }
            if (y < y_min || y > y_max) {
                throw std::runtime_error(
                    "Interpolation2D::operator():y = " + std::to_string(y) +
                    " is outside interpolation area (" + std::to_string(y_min) +
                    ", " + std::to_string(y_max) + ").");
            }
            return gsl_spline2d_eval(spline, x, y, x_acc, y_acc);
#endif
        }

        ~Interpolation2D() {
            gsl_spline2d_free(spline);
            gsl_interp_accel_free(x_acc);
            gsl_interp_accel_free(y_acc);
        }
};

#endif /* ifndef INTERPOLATION_HPP */
