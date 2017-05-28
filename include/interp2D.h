#pragma once

#include <vector>
#include <array>
#include <string>
#include <boost/multi_array.hpp>
#include <netcdf>

using boost::multi_array;
using boost::extents;
using namespace netCDF;

struct axis {
    uint32_t nx;
    double   dx;
    std::vector<double> x;

    void reset() {
        nx = 0;
        dx = 0;
        x.clear();
    }
};

typedef std::array<axis, 2>           axes2D;
typedef boost::multi_array<double, 2> data2D;

class netcdf_reader {

    NcFile ncdf_file;
    std::vector<std::string> ncdf_dims, ncdf_vars;

    axes2D a2d;
    std::vector<data2D> d2d;

    void                    load_dims();
    void                    load_data();

  public:
    netcdf_reader(const std::string &nc_file,
                  const std::vector<std::string> &dims,
                  const std::vector<std::string> &vars);

    virtual ~netcdf_reader() {}


    void                    load();

    axes2D                  get_dims();
    std::vector<data2D>     get_data();
};

class interp2D {
  public:
    axes2D                  axes;
    std::vector<data2D>     data;

  public:

    interp2D() {};
    interp2D(const std::string &ncdf_file);
    interp2D(const axes2D &input_axis, const std::vector<data2D> &input_data);

    virtual ~interp2D() {}

    void load_from_netcdf(const std::string &nc_file);

    void set_data(const axes2D &input_axis, const std::vector<data2D> &input_data);

    bool   check_point(double xval, double yval);
    double interpolate(int col_id, double xval, double yval);

};

