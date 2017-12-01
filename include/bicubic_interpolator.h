#pragma once

#include <vector>
#include <array>
#include <string>
#include <boost/multi_array.hpp>
#include <netcdf>

#include "axis.h"
#include "netcdf_reader.h"

using boost::multi_array;
using boost::extents;

#define CLAMP(v, min, max) if (v < min) { v = min; } else if (v > max) { v = max; }

typedef boost::multi_array<double, 4> interp_matrix_type;

class bicubic_interpolator
{
private:
    axes2D                  a;
    std::vector<data2D>     f;

    std::vector<interp_matrix_type> c;

    double get_value_clamped(size_t id, size_t yidx, size_t xidx);
    uint8_t check_point ( double x, double y );
    void make_c();
    std::array<double, 4> mmul(const std::array<double, 4>& x);

    double interp_val(size_t idx, double x, double y);

public:

    bicubic_interpolator () {}
    bicubic_interpolator ( const axes2D &input_axis, const std::vector<data2D> &input_data );
    bicubic_interpolator ( const std::string &ncdf_file );

    virtual ~bicubic_interpolator() {}

    double operator () (size_t idx, double x, double y);
    void  set_data ( const axes2D &input_axis, const std::vector<data2D> &input_data );
    void load_from_netcdf ( const std::string &nc_file );

};

