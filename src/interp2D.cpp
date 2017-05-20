#include <bits/stdc++.h>
#include <boost/multi_array.hpp>
#include <plog/Log.h>
#include <netcdf>

#include "interp2D.h"

using boost::multi_array;
using boost::extents;
using namespace netCDF;

netcdf_reader::netcdf_reader(const std::string& nc_file, 
                const std::vector<std::string>& dims, 
                const std::vector<std::string>& vars) : 
  ncdf_file(nc_file, NcFile::read),
  ncdf_dims(dims), ncdf_vars(vars)
{

}


void
netcdf_reader::load_dims()
{
  std::multimap<std::string, NcDim> dims = ncdf_file.getDims();
  
  for(auto &d : dims)
  {
    auto didx = -1;
    for(auto i = 0; i < ncdf_dims.size(); i++)
    {
      if (d.first == ncdf_dims[i]) didx = i;
    }

    if (didx == -1)
    {
      std::cerr<<"No dim named "<< d.first << std::endl;
      exit(1);
    }

    a2d[didx].nx = d.second.getSize(); 

    a2d[didx].x.resize(a2d[didx].nx);

    auto buf = std::make_unique<double[]>(a2d[didx].nx);
    auto var = ncdf_file.getVar(d.first.c_str());

    if (var.isNull()) 
    {
      std::cerr << "FAILED TO GET "<< d.first << std::endl;
      exit(1);
    }
    else
    {
      auto pbuf = (double*)buf.get();
      var.getVar(pbuf);
      a2d[didx].x.assign(pbuf, pbuf + a2d[didx].nx);
    }
    a2d[didx].dx = (a2d[didx].x.back() - a2d[didx].x.front()) / ((double) a2d[didx].nx);
  }
}


void
netcdf_reader::load_data()
{
  d2d.resize(ncdf_vars.size()); 
 
  for(auto &rd : d2d)
    rd.resize(extents[a2d[0].nx][a2d[1].nx]);

  auto buflen = a2d[0].nx * a2d[1].nx;

  auto buffer = std::make_unique<double[]>(buflen);
  auto pbuffer = (double*)buffer.get();

  for(auto i = 0 ; i < ncdf_vars.size(); i++)
  {
    NcVar ncdf_data = ncdf_file.getVar(ncdf_vars[i]);

    ncdf_data.getVar(pbuffer); 
    d2d[i].assign(pbuffer, pbuffer + buflen);
  }
} 

axes2D
netcdf_reader::get_dims() { return a2d; }

std::vector<data2D>
netcdf_reader::get_data() { return d2d; }

interp2D::interp2D(const std::string& ncdf_file)
{
  load_from_netcdf(ncdf_file);
}

interp2D::interp2D(const axes2D& input_axis, const std::vector<data2D>& input_data)
{
  set_data(input_axis, input_data);
}

void
interp2D::set_data(const axes2D& input_axis, const std::vector<data2D>& input_data)
{
  axes = input_axis;
  data.assign(input_data.begin(),input_data.end());
}
void
interp2D::load_from_netcdf(const std::string& nc_file)
{
  std::vector<std::string> dims = {"temperature", "saturation"};
  std::vector<std::string> cols = {"nucleation_rate", "critical_size"};

  netcdf_reader nr(nc_file, dims, cols);

  set_data(nr.get_dims(), nr.get_data());  
}

double
interp2D::interpolate(int col_id, double xval, double yval)
{

  auto cr = [=] (auto val, auto dim) 
    { 
      if (val < axes[dim].x.front() || val >= axes[dim].x.back())
        return false;
      return true;
    };
  
  auto get_idx = [=] (auto val, auto dim)
    {
      return (int)((val - axes[dim].x[dim]) / (axes[dim].dx));
    };

  if (!cr(xval, 0) || !cr(yval, 1)) 
  {
    LOGE << "given value is outside of interpolate axis";
    return 0.0;      
  }

  int ix = get_idx(xval, 0);
  int iy = get_idx(yval, 1);

//  std::cout << axes[0].dx <<" , "<<axes[1].dx<<std::endl;

  double xp = (xval - axes[0].x[ix]) / axes[0].dx;
  double yp = (yval - axes[1].x[iy]) / axes[1].dx;
 
  double f0 = data[col_id][ix][iy];
  double f1 = data[col_id][ix + 1][iy];
  double f2 = data[col_id][ix][iy + 1];
  double f3 = data[col_id][ix + 1][iy + 1];

  double r1 = f0 * (1.0 - xp) * (1.0 - yp)
              + f1 * xp * (1.0 - yp)
              + f2 * (1.0 - xp) * yp
              + f3 * xp * yp;

  return r1;
}

