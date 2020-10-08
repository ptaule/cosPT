/*
   interpolation.cpp

   Created by Petter Taule on 04.09.2020
   Copyright (c) 2020 Petter Taule. All rights reserved.
*/

#include <vector>
#include <stdexcept>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include "../include/io.hpp"
#include "../include/interpolation.hpp"

void Interpolation1D::initialize(const Vec1D<double>& x, const Vec1D<double>& y)
{
    if (x.size() != y.size()) {
        throw(std::invalid_argument(
            "Dimensions of x and y vectors are not equal."));
    }

    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(type, x.size());
    gsl_spline_init(spline, x.data(), y.data(), x.size());
}



Interpolation1D::Interpolation1D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const gsl_interp_type* type
        ) : type(type)
{
    initialize(x,y);
}

Interpolation1D::Interpolation1D(
        const Vec1D<double>& x,
        const Vec1D<double>& y
        ) : Interpolation1D(x, y, gsl_interp_cspline)
{}



Interpolation1D::Interpolation1D(
        const std::string& filename,
        const gsl_interp_type* type
        ) : type(type)
{
    Vec2D<double> columns;
    read_columns_from_file(filename, 2, columns);
    initialize(columns.at(0), columns.at(1));
}



Interpolation1D::Interpolation1D(const std::string& filename)
    : Interpolation1D(filename, gsl_interp_cspline)
{}



Interpolation1D::Interpolation1D(Interpolation1D&& other) {
    spline = other.spline;
    acc    = other.acc;
    type   = other.type;

    other.spline = nullptr;
    other.acc    = nullptr;
    other.type   = nullptr;
}



Interpolation1D& Interpolation1D::operator=(Interpolation1D&& other)
{
    if (this != &other) {
        spline = other.spline;
        acc = other.acc;
        type = other.type;

        other.spline = nullptr;
        other.acc    = nullptr;
        other.type   = nullptr;
    }
    return *this;
}



void Interpolation2D::initialize(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const Vec1D<double>& z
        ) 
{
    if (x.size() * y.size() != z.size()) {
        throw(std::invalid_argument(
            "Length of z vector does not equal x.size() * y.size()."));
    }

    x_acc = gsl_interp_accel_alloc();
    y_acc = gsl_interp_accel_alloc();
    spline = gsl_spline2d_alloc(type, x.size(), y.size());
    gsl_spline2d_init(spline, x.data(), y.data(), z.data(), x.size(), y.size());

}



Interpolation2D::Interpolation2D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const Vec1D<double>& z,
        const gsl_interp2d_type* type
        ) : type(type)
{
    initialize(x, y, z);
}



Interpolation2D::Interpolation2D(
        const Vec1D<double>& x,
        const Vec1D<double>& y,
        const Vec1D<double>& z
        ) : Interpolation2D(x, y, z, gsl_interp2d_bicubic)
{}



Interpolation2D::Interpolation2D(
        const std::string& x_grid_file,
        const std::string& y_grid_file,
        const std::string& data_file,
        const gsl_interp2d_type* type
        ) : type(type)
{
    Vec2D<double> x;
    Vec2D<double> y;
    Vec1D<double> z;
    read_columns_from_file(x_grid_file, 1, x);
    read_columns_from_file(y_grid_file, 1, y);

    read_data_grid_from_file(data_file, z, x.at(0).size(), y.at(0).size());

    initialize(x.at(0), y.at(0), z);
}



Interpolation2D::Interpolation2D(
        const std::string& x_grid_file,
        const std::string& y_grid_file,
        const std::string& data_file
        ) : Interpolation2D(x_grid_file, y_grid_file, data_file,
            gsl_interp2d_bicubic)
{}



Interpolation2D::Interpolation2D(Interpolation2D&& other)
{
    x_acc = other.x_acc;
    y_acc = other.y_acc;
    spline = other.spline;
    type = other.type;

    other.x_acc = nullptr;
    other.y_acc = nullptr;
    other.spline = nullptr;
    other.type = nullptr;
}



Interpolation2D& Interpolation2D::operator=(Interpolation2D&& other)
{
    if (this != &other) {
        spline = other.spline;
        x_acc = other.x_acc;
        y_acc = other.y_acc;
        type = other.type;

        other.spline = nullptr;
        other.x_acc  = nullptr;
        other.y_acc  = nullptr;
        other.type   = nullptr;
    }
    return *this;
}
