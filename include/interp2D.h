#pragma once

#include <bits/stdc++.h>
#include <boost/multi_array.hpp>
#include <netcdf>

using boost::multi_array;
using boost::extents;
using namespace netCDF;

struct axis
{
  uint32_t nx;
  double   dx;
  std::vector<double> x;

  void reset() { nx = 0; dx = 0; x.clear(); }
};

typedef std::array<axis, 2>           axes2D;
typedef boost::multi_array<double, 2> data2D;

class netcdf_reader
{

  NcFile ncdf_file;
  std::vector<std::string> ncdf_dims, ncdf_vars;
  
public:
  netcdf_reader(const std::string& nc_file, 
                const std::vector<std::string>& dims, 
                const std::vector<std::string>& vars);

  virtual ~netcdf_reader() {}

  
  axes2D                  get_dims();
  std::array<data2D, 2>   get_data(const axes2D& ax);
};

class interp2D
{
  axes2D                  axes;
  std::vector<data2D>     data;

public:

  interp2D(const axes2D& input_axis, const std::vector<data2D>& data);

  virtual ~interp2D() {}

  double interpolate(int col_id, double xval, double yval);

};

