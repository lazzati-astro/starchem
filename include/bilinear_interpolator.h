#pragma once

#include <vector>
#include <array>
#include <string>
#include <boost/multi_array.hpp>

#include "axis.h"
#include "netcdf_reader.h"

using boost::multi_array;
using boost::extents;

class bilinear_interpolator
{
public:
    axes2D                  axes;
    std::vector<data2D>     data;

public:

    bilinear_interpolator() {};
    bilinear_interpolator ( const std::string &ncdf_file );
    bilinear_interpolator ( const axes2D &input_axis, const std::vector<data2D> &input_data );

    virtual ~bilinear_interpolator() {}

    void load_from_netcdf ( const std::string &nc_file );

    void set_data ( const axes2D &input_axis, const std::vector<data2D> &input_data );

    bool   check_point ( double xval, double yval );
    double operator () ( int col_id, double xval, double yval );

};

