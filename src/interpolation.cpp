/*
   interpolation.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <vector>
#include <algorithm>
#include <stdexcept>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "../include/io.hpp"
#include "../include/interpolation.hpp"

void Interpolation1D::initialize(
        Vec1D<double> x,
        Vec1D<double> y,
        double factor,
        bool padding_zeros
        )
{
    if (x.size() != y.size()) {
        throw(std::invalid_argument(
            "Dimensions of x and y vectors are not equal."));
    }

    if (padding_zeros) {
        /* y(x = 0) = 0 */
        x.insert(x.begin(), 0);
        y.insert(y.begin(), 0);

        /* y(10 * x_max) = 0 */
        x.push_back(10 * x.back());
        y.push_back(0);
    }

    /* If factor != 1, multiply y-values by factor */
    if (factor != 1) {
        std::transform(y.begin(), y.end(), y.begin(),
                [factor](double& d) -> double {return factor * d;});
    }

    x_min = x.front();
    x_max = x.back();

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(type, x.size());
    gsl_spline_init(spline, x.data(), y.data(), x.size());
}



Interpolation1D::Interpolation1D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        double factor,
        bool padding_zeros,
        const gsl_interp_type* type
        ) : type(type)
{
    initialize(x, y, factor, padding_zeros);
}



Interpolation1D::Interpolation1D(
        const std::string& filename,
        double factor,
        bool padding_zeros,
        const gsl_interp_type* type
        ) : type(type)
{
    Vec2D<double> columns;
    read_columns_from_file(filename, 2, columns);
    initialize(columns.at(0), columns.at(1), factor, padding_zeros);
}



Interpolation1D::Interpolation1D(Interpolation1D&& other)
    : x_min (std::exchange(other.x_min, 0)),
    x_max(std::exchange(other.x_max, 0)),
    spline(std::exchange(other.spline, nullptr)),
    acc(std::exchange(other.acc, nullptr)),
    type(std::exchange(other.type, nullptr))
{}



Interpolation1D& Interpolation1D::operator=(Interpolation1D&& other)
{
    if (this != &other) {
        x_min  = std::exchange(other.x_min, 0);
        x_max  = std::exchange(other.x_max, 0);
        acc    = std::exchange(other.acc, nullptr);
        spline = std::exchange(other.spline, nullptr);
        type   = std::exchange(other.type, nullptr);
    }
    return *this;
}



void Interpolation2D::initialize(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        Vec1D<double> z,
        double factor
        )
{
    if (x.size() * y.size() != z.size()) {
        throw(std::invalid_argument(
            "Length of z vector does not equal x.size() * y.size()."));
    }

    /* If factor != 1, multiply y-values by factor */
    if (factor != 1) {
        std::transform(z.begin(), z.end(), z.begin(),
                [factor](double& d) -> double {return factor * d;});
    }

    x_min = x.front();
    x_max = x.back();
    y_min = y.front();
    y_max = y.back();

    x_acc = gsl_interp_accel_alloc();
    y_acc = gsl_interp_accel_alloc();
    spline = gsl_spline2d_alloc(type, x.size(), y.size());
    gsl_spline2d_init(spline, x.data(), y.data(), z.data(), x.size(), y.size());

}



Interpolation2D::Interpolation2D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const Vec1D<double>& z,
        double factor,
        const gsl_interp2d_type* type
        ) : type(type)
{
    initialize(x, y, z, factor);
}



Interpolation2D::Interpolation2D(
        const std::string& x_grid_file,
        const std::string& y_grid_file,
        const std::string& data_file,
        double factor,
        const gsl_interp2d_type* type
        ) : type(type)
{
    Vec2D<double> x;
    Vec2D<double> y;
    Vec1D<double> z;
    read_columns_from_file(x_grid_file, 1, x);
    read_columns_from_file(y_grid_file, 1, y);

    read_data_grid_from_file(data_file, z, x.at(0).size(), y.at(0).size());

    initialize(x.at(0), y.at(0), z, factor);
}



Interpolation2D::Interpolation2D(Interpolation2D&& other) noexcept
    : x_min (std::exchange(other.x_min, 0)),
    x_max(std::exchange(other.x_max, 0)),
    y_min(std::exchange(other.y_min, 0)),
    y_max(std::exchange(other.y_max, 0)),
    spline(std::exchange(other.spline, nullptr)),
    x_acc(std::exchange(other.x_acc, nullptr)),
    y_acc(std::exchange(other.y_acc, nullptr)),
    type(std::exchange(other.type, nullptr))
{}



Interpolation2D& Interpolation2D::operator=(Interpolation2D&& other)
{
    if (this != &other) {
        x_min  = std::exchange(other.x_min, 0);
        x_max  = std::exchange(other.x_max, 0);
        y_min  = std::exchange(other.y_min, 0);
        y_max  = std::exchange(other.y_max, 0);

        x_acc  = std::exchange(other.x_acc, nullptr);
        y_acc  = std::exchange(other.y_acc, nullptr);
        spline = std::exchange(other.spline, nullptr);
        type   = std::exchange(other.type, nullptr);
    }
    return *this;
}
