#pragma once

#include <vector>
#include <string>
#include <netcdf>

#include "axis.h"

using namespace netCDF;
using boost::extents;

class netcdf_reader
{

    NcFile ncdf_file;
    std::vector<std::string> ncdf_dims, ncdf_vars;

    axes2D a2d;
    std::vector<data2D> d2d;

    void load_dims();
    void load_data();

public:

    netcdf_reader ( const std::string &nc_file,
                    const std::vector<std::string> &dims,
                    const std::vector<std::string> &vars ) ;

    axes2D get_dims();

    std::vector<data2D> get_data();
};


