#include <vector>
#include <map>
#include <string>
#include <netcdf>
#include <boost/multi_array.hpp>

#include "axis.h"
#include "netcdf_reader.h"

using namespace netCDF;
using boost::extents;

void
netcdf_reader::load_dims()
{
    std::multimap<std::string, NcDim> dims = ncdf_file.getDims();

    for ( auto &d : dims )
    {
        auto didx = -1;
        for ( auto i = 0; i < ncdf_dims.size(); i++ )
        {
            if ( d.first == ncdf_dims[i] ) didx = i;
        }

        if ( didx == -1 )
        {
            std::cerr << "No dim named " << d.first << std::endl;
            exit ( 1 );
        }

        a2d[didx].nx = d.second.getSize();

        a2d[didx].x.resize ( a2d[didx].nx );

        auto buf = std::make_unique<double[]> ( a2d[didx].nx );
        auto var = ncdf_file.getVar ( d.first.c_str() );

        if ( var.isNull() )
        {
            std::cerr << "FAILED TO GET " << d.first << std::endl;
            exit ( 1 );
        }
        else
        {
            auto pbuf = ( double * ) buf.get();
            var.getVar ( pbuf );
            a2d[didx].x.assign ( pbuf, pbuf + a2d[didx].nx );
        }
        a2d[didx].dx = ( a2d[didx].x.back() - a2d[didx].x.front() ) / ( ( double ) ( a2d[didx].nx - 1 ) );
        a2d[didx].idx = 1. / a2d[didx].dx;
    }
}


void netcdf_reader::load_data()
{
    d2d.resize ( ncdf_vars.size() );

    for ( auto &rd : d2d )
        rd.resize ( extents[a2d[0].nx][a2d[1].nx] );

    auto buflen = a2d[0].nx * a2d[1].nx;

    for ( auto i = 0 ; i < ncdf_vars.size(); i++ )
    {
        auto buffer = std::make_unique<double[]> ( buflen );
        auto pbuffer = ( double * ) buffer.get();

        NcVar ncdf_data = ncdf_file.getVar ( ncdf_vars[i] );

        ncdf_data.getVar ( pbuffer );
        d2d[i].assign ( pbuffer, pbuffer + buflen );
    }
}


netcdf_reader::netcdf_reader ( const std::string &nc_file,
                               const std::vector<std::string> &dims,
                               const std::vector<std::string> &vars ) :
    ncdf_file ( nc_file, NcFile::read ),
    ncdf_dims ( dims ), ncdf_vars ( vars )
{
    load_dims();
    load_data();
}


axes2D netcdf_reader::get_dims()
{
    return a2d;
}

std::vector<data2D> netcdf_reader::get_data()
{
    return d2d;
}


