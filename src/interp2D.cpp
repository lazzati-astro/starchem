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


axes2D
netcdf_reader::get_dims()
{
  axes2D axes;

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

    axes[didx].nx = d.second.getSize(); 

    axes[didx].x.resize(axes[didx].nx);

    auto buf = std::make_unique<double[]>(axes[didx].nx);
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
      axes[didx].x.assign(pbuf, pbuf + axes[didx].nx);
    }
    axes[didx].dx = (axes[didx].x.back() - axes[didx].x.front()) / ((double) axes[didx].nx);
  }
  return axes;
}


std::array<data2D, 2>
netcdf_reader::get_data(const axes2D& ax)
{
  std::array<data2D, 2> rdata;
  
  for(auto &rd : rdata)
    rd.resize(extents[ax[0].nx][ax[1].nx]);

  auto buflen = ax[0].nx * ax[1].nx;

  auto buffer = std::make_unique<double[]>(buflen);
  auto pbuffer = (double*)buffer.get();

  for(auto i = 0 ; i < ncdf_vars.size(); i++)
  {
    NcVar ncdf_data = ncdf_file.getVar(ncdf_vars[i]);

    ncdf_data.getVar(pbuffer); 
    rdata[i].assign(pbuffer, pbuffer + buflen);
  }
  return rdata;
} 

interp2D::interp2D(const axes2D& input_axis, const std::vector<data2D>& input_data)
{
  axes = input_axis;
  data.assign(input_data.begin(),input_data.end());
}

double
interp2D::interpolate(int col_id, double xval, double yval)
{
  int ix = (int)((xval - axes[0].x[0]) / (axes[0].dx));
  int iy = (int)((yval - axes[1].x[0]) / (axes[1].dx));

  std::cout << axes[0].dx <<" , "<<axes[1].dx<<std::endl;

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

