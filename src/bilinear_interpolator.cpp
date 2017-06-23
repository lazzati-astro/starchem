#include <vector>
#include <string>
#include <map>

#include <boost/multi_array.hpp>
#include <plog/Log.h>
#include <netcdf>

#include "axis.h"
#include "netcdf_reader.h"
#include "bilinear_interpolator.h"

using boost::multi_array;
using boost::extents;
using namespace netCDF;

bilinear_interpolator::bilinear_interpolator ( const std::string &ncdf_file )
{
    load_from_netcdf ( ncdf_file );
}

bilinear_interpolator::bilinear_interpolator ( const axes2D &input_axis, const std::vector<data2D> &input_data )
{
    std::cout << input_axis.size() << std::endl;
    set_data ( input_axis, input_data );
}

void
bilinear_interpolator::set_data ( const axes2D &input_axis, const std::vector<data2D> &input_data )
{
    axes = input_axis;
    data.assign ( input_data.begin(), input_data.end() );
}
void
bilinear_interpolator::load_from_netcdf ( const std::string &nc_file )
{
    std::vector<std::string> dims = {"temperature", "saturation"};
    std::vector<std::string> cols = {"nucleation_rate", "critical_size"};

    netcdf_reader nr ( nc_file, dims, cols );

    set_data ( nr.get_dims(), nr.get_data() );
}

bool
bilinear_interpolator::check_point ( double xval, double yval )
{
    if ( xval <= axes[0].x.front() ||
            xval >= axes[0].x.back()  ||
            yval <= axes[1].x.front() ||
            yval >= axes[1].x.back() ) return false;
    return true;
}

double
bilinear_interpolator::operator () ( int col_id, double xval, double yval )
{

    if ( !check_point ( xval, yval ) )
    {
        return -1.0;
    }

    auto const itx = std::lower_bound ( axes[0].x.begin(), axes[0].x.end(), xval ) - 1;
    auto const ity = std::lower_bound ( axes[1].x.begin(), axes[1].x.end(), yval ) - 1;

    auto const ix = std::distance ( axes[0].x.begin(), itx );
    auto const iy = std::distance ( axes[1].x.begin(), ity );

    double xp = ( xval - *itx ) * axes[0].idx;
    double yp = ( yval - *ity ) * axes[1].idx;

    auto ixp1 = ix + 1;
    ixp1 = ( ixp1 >= axes[0].x.size() ? ix : ixp1 );

    auto iyp1 = iy + 1;
    iyp1 = ( iyp1 >= axes[1].x.size() ? iy : iyp1 );

    if (iy < 0)
    {
        std::cout << "yval = " << yval << std::endl;
    }
    assert (ix >= 0);
    assert (iy >= 0);
    assert (ixp1 >= 0);
    assert (iyp1 >= 0);


    double f0 = data[col_id][iy][ix];
    double f1 = data[col_id][iyp1][ix];
    double f2 = data[col_id][iy][ixp1];
    double f3 = data[col_id][iyp1][ixp1];

    double r1 = f0 * ( 1.0 - xp ) * ( 1.0 - yp )
                + f1 * xp * ( 1.0 - yp )
                + f2 * ( 1.0 - xp ) * yp
                + f3 * xp * yp;

    return r1;
}

