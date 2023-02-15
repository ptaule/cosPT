/*
   interpolation.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "../include/io.hpp"
#include "../include/interpolation.hpp"

#if (__cplusplus < 201402L)
template <typename T, typename U>
T exchange(T &storage, U &&value) {
    T out = std::move(storage);
    storage = (T)std::move(value);
    return out;
}
#else
using std::exchange;
#endif

using std::size_t;

template <class T>
using Vec1D = std::vector<T>;
template <class T>
using Vec2D = std::vector<std::vector<T>>;

void Interpolation1D::ensure_increasing_x(
        std::vector<double>& x,
        std::vector<double>& y
        )
{
    /* Are the x-values increasing/decreasing? */
    bool x_increasing = true;
    if (x.at(1) < x.at(0)) {
        x_increasing = false;
    }

    /* Check that the rest of the elements are also increasing/decreasing */
    for (size_t i = 2; i < x.size(); ++i) {
        if (x_increasing && (x.at(i) <= x.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation1D::ensure_increasing_x(): x-values are "
                "not strictly increasing or decreasing.");
        }
        if (!x_increasing && (x.at(i) >= x.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation1D::ensure_increasing_x(): x-values are "
                "not strictly increasing or decreasing.");
        }
    }

    /* Reverse x- and y-values if x_increasing is false (the GSL interpolation
     * code requires increasing x values). */
    if (!x_increasing) {
        std::reverse(x.begin(), x.end());
        std::reverse(y.begin(), y.end());
    }
}



void Interpolation1D::initialize(
        Vec1D<double> x,
        Vec1D<double> y,
        double factor
        )
{
    if (x.size() != y.size()) {
        throw(std::invalid_argument("Interpolation1D::initialize(): Dimensions "
                                    "of x and y vectors are not equal."));
    }

    ensure_increasing_x(x,y);
    x_min = x.front();
    x_max = x.back();

    /* If factor != 1, multiply y-values by factor */
    if (factor != 1) {
        std::transform(y.begin(), y.end(), y.begin(),
                [factor](double& d) -> double {return factor * d;});
    }

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(type, x.size());
    gsl_spline_init(spline, x.data(), y.data(), x.size());
}



Interpolation1D::Interpolation1D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        double factor,
        const gsl_interp_type* type
        ) : type(type)
{
    initialize(x, y, factor);
}



Interpolation1D::Interpolation1D(
        const std::string& filename,
        double factor,
        const gsl_interp_type* type
        ) : type(type)
{
    Vec2D<double> data;
    read_delimited_file(filename, data);

    /* data is row-by-row, want to convert it to two columns x,y */
    Vec1D<double> x(data.size());
    Vec1D<double> y(data.size());

    for (size_t i = 0; i < data.size(); ++i) {
        try {
            x.at(i) = data.at(i).at(0);
            y.at(i) = data.at(i).at(1);
        }
        catch (std::out_of_range& e) {
            std::cerr
                << "Interpolation1D::Interpolation1D(): Encountered "
                   "out_of_range exception after reading \""
                << filename
                << "\" for interpolation. This file should have two columns."
                << std::endl;
            throw e;
        }
    }
    initialize(x, y, factor);
}



Interpolation1D::Interpolation1D(Interpolation1D&& other)
    : x_min (exchange(other.x_min, 0)),
    x_max(exchange(other.x_max, 0)),
    spline(exchange(other.spline, nullptr)),
    acc(exchange(other.acc, nullptr)),
    type(exchange(other.type, nullptr))
{}



Interpolation1D& Interpolation1D::operator=(Interpolation1D&& other)
{
    if (this != &other) {
        x_min  = exchange(other.x_min, 0);
        x_max  = exchange(other.x_max, 0);
        acc    = exchange(other.acc, nullptr);
        spline = exchange(other.spline, nullptr);
        type   = exchange(other.type, nullptr);
    }
    return *this;
}



void Interpolation2D::ensure_increasing_x_y(
        std::vector<double>& x,
        std::vector<double>& y,
        std::vector<std::vector<double>>& z
        )
{
    /* Are the x-values increasing/decreasing? */
    bool x_increasing = true;
    if (x.at(1) < x.at(0)) {
        x_increasing = false;
    }

    /* Check that the rest of the elements are also increasing/decreasing */
    for (size_t i = 2; i < x.size(); ++i) {
        if (x_increasing && (x.at(i) <= x.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation2D::ensure_increasing_x_y(): x-values are "
                "not strictly increasing or decreasing.");
        }
        if (!x_increasing && (x.at(i) >= x.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation2D::ensure_increasing_x_y(): x-values are "
                "not strictly increasing or decreasing.");
        }
    }

    /* Reverse x- and z-values (first dimension) if x_increasing is false (the
     * GSL interpolation code requires increasing x values). */
    if (!x_increasing) {
        std::reverse(x.begin(), x.end());
        std::reverse(z.begin(), z.end());
    }

    /* Are the y-values increasing/decreasing? */
    bool y_increasing = true;
    if (y.at(1) < y.at(0)) {
        y_increasing = false;
    }

    /* Check that the rest of the elements are also increasing/decreasing */
    for (size_t i = 2; i < y.size(); ++i) {
        if (y_increasing && (y.at(i) <= y.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation2D::ensure_increasing_x_y(): y-values are "
                "not strictly increasing or decreasing.");
        }
        if (!y_increasing && (y.at(i) >= y.at(i-1))) {
            throw std::invalid_argument(
                "Interpolation2D::ensure_increasing_x_y(): y-values are "
                "not strictly increasing or decreasing.");
        }
    }

    /* Reverse y- and z-values (second dimension) if y_increasing is false
     * (the GSL interpolation code requires increasing y values). */
    if (!y_increasing) {
        std::reverse(y.begin(), y.end());
        for (auto& el : z) {
            std::reverse(el.begin(), el.end());
        }
    }
}



void Interpolation2D::initialize(
        Vec1D<double> x,
        Vec1D<double> y,
        Vec2D<double> z,
        double factor
        )
{
    /* Check that dimensions matches */
    if (x.size() != z.size()) {
        throw(std::invalid_argument(
            "Interpolation2D::initialize(): Dimension mismatch between x-vector "
            "and first dimension of z-matrix."));
    }
    for (auto& el : z) {
        if (y.size() != el.size()) {
            throw(std::invalid_argument("Interpolation2D::initialize(): "
                                        "Dimension mismatch between y-vector "
                                        "and second dimension of z-matrix."));
        }
    }

    ensure_increasing_x_y(x, y, z);
    x_min = x.front();
    x_max = x.back();
    y_min = y.front();
    y_max = y.back();

    /* If factor != 1, multiply z-values by factor */
    if (factor != 1) {
        for (auto& el : z) {
            std::transform(el.begin(), el.end(), el.begin(),
                    [factor](double& d) -> double {return factor * d;});
        }
    }

    Vec1D<double> z_1d(x.size() * y.size());

    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t j = 0; j < y.size(); ++j) {
            z_1d.at(idx2d_to_idx1d(i,j, x.size())) = z.at(i).at(j);
        }
    }

    x_acc = gsl_interp_accel_alloc();
    y_acc = gsl_interp_accel_alloc();
    spline = gsl_spline2d_alloc(type, x.size(), y.size());
    gsl_spline2d_init(spline, x.data(), y.data(), z_1d.data(), x.size(), y.size());
}



Interpolation2D::Interpolation2D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const Vec2D<double>& z,
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
    Vec2D<double> x_2d;
    Vec2D<double> y_2d;
    Vec2D<double> z;
    read_delimited_file(x_grid_file, x_2d);
    read_delimited_file(y_grid_file, y_2d);
    read_delimited_file(data_file,   z);

    /* x_2d and y_2d is row-by-row, want to convert each to a column */
    Vec1D<double> x(x_2d.size());
    Vec1D<double> y(y_2d.size());

    for (size_t i = 0; i < x_2d.size(); ++i) {
        if (x_2d.at(i).size() != 1) {
            std::cerr << "Interpolation2D::Interpolation2D(): Warning: Found "
                         "more than one column in \""
                      << x_grid_file
                      << "\", only first column used for interpolation."
                      << std::endl;
        }
        x.at(i) = x_2d.at(i).at(0);
    }
    for (size_t i = 0; i < y_2d.size(); ++i) {
        if (y_2d.at(i).size() != 1) {
            std::cerr << "Interpolation2D::Interpolation2D(): Warning: Found "
                         "more than one column in \""
                      << y_grid_file
                      << "\", only first column used for interpolation."
                      << std::endl;
        }
        y.at(i) = y_2d.at(i).at(0);
    }
    initialize(x, y, z, factor);
}



Interpolation2D::Interpolation2D(Interpolation2D&& other) noexcept
    : x_min (exchange(other.x_min, 0)),
    x_max(exchange(other.x_max, 0)),
    y_min(exchange(other.y_min, 0)),
    y_max(exchange(other.y_max, 0)),
    spline(exchange(other.spline, nullptr)),
    x_acc(exchange(other.x_acc, nullptr)),
    y_acc(exchange(other.y_acc, nullptr)),
    type(exchange(other.type, nullptr))
{}



Interpolation2D& Interpolation2D::operator=(Interpolation2D&& other)
{
    if (this != &other) {
        x_min  = exchange(other.x_min, 0);
        x_max  = exchange(other.x_max, 0);
        y_min  = exchange(other.y_min, 0);
        y_max  = exchange(other.y_max, 0);

        x_acc  = exchange(other.x_acc, nullptr);
        y_acc  = exchange(other.y_acc, nullptr);
        spline = exchange(other.spline, nullptr);
        type   = exchange(other.type, nullptr);
    }
    return *this;
}
