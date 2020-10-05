/*
   interpolation.hpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
   */

#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include <vector>
#include <string>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "utilities.hpp"

/* Wrapper classes around GSL 1d and 2d spline
 * interpolation */

class Interpolation1D {
    private:
        gsl_spline* spline;
        gsl_interp_accel* acc;
        const gsl_interp_type* type;

        void initialize(const Vec1D<double>& x, const Vec1D<double>& y);
    public:
        Interpolation1D() : spline(nullptr), acc(nullptr), type(nullptr) {}
        Interpolation1D(
                const Vec1D<double>& x,
                const Vec1D<double>& y,
                const gsl_interp_type* type
                );
        Interpolation1D(
                const Vec1D<double>& x,
                const Vec1D<double>& y
                );
        Interpolation1D(const std::string& filename, const gsl_interp_type* type);
        Interpolation1D(const std::string& filename);
        Interpolation1D(const Interpolation1D&) = delete;
        Interpolation1D(Interpolation1D&& other);


        double eval(double x) const {
            return gsl_spline_eval(spline, x, acc);
        }

        ~Interpolation1D() {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        }
};


class Interpolation2D {
    private:
        gsl_spline2d* spline;
        gsl_interp_accel* x_acc;
        gsl_interp_accel* y_acc;
        const gsl_interp2d_type* type;

        void initialize(
                const Vec1D<double>& x,
                const Vec1D<double>& y,
                const Vec1D<double>& z
                );
    public:
        Interpolation2D() :
            spline(nullptr), x_acc(nullptr), y_acc(nullptr), type(nullptr) {}
        Interpolation2D(
                const Vec1D<double>& x,
                const Vec1D<double>& y,
                const Vec1D<double>& z,
                const gsl_interp2d_type* type
                );
        Interpolation2D(
                const Vec1D<double>& x,
                const Vec1D<double>& y,
                const Vec1D<double>& z
                );
        Interpolation2D(
                const std::string& x_grid_file,
                const std::string& y_grid_file,
                const std::string& data_file,
                const gsl_interp2d_type* type
                );
        Interpolation2D(
                const std::string& x_grid_file,
                const std::string& y_grid_file,
                const std::string& data_file
                );
        Interpolation2D(const Interpolation2D&) = delete;
        Interpolation2D(Interpolation2D&&);


        double eval(double x, double y) const {
            return gsl_spline2d_eval(spline, x, y, x_acc, y_acc);
        }

        ~Interpolation2D() {
            gsl_spline2d_free(spline);
            gsl_interp_accel_free(x_acc);
            gsl_interp_accel_free(y_acc);
        }
};

#endif /* ifndef INTERPOLATION_HPP */
