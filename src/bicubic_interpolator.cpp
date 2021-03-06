#include <vector>
#include <array>
#include <string>
#include <boost/multi_array.hpp>
#include <netcdf>

#include "axis.h"
#include "netcdf_reader.h"
#include "bicubic_interpolator.h"

using boost::multi_array;
using boost::extents;

double
bicubic_interpolator::get_value_clamped( size_t id, size_t yidx, size_t xidx )
{
    CLAMP(xidx, 0, a[0].nx - 1);
    CLAMP(yidx, 0, a[1].nx - 1);
    return f[id][yidx][xidx];
}

uint8_t
bicubic_interpolator::check_point ( double x, double y )
{
    if ( x <= a[0].x.front() ) return false;
    if ( x > a[0].x.back() ) return false;
    if ( y <= a[1].x.front()) return false;
    if ( y > a[1].x.back() ) return false;
    return true;
}

void
bicubic_interpolator::make_c()
{

    std::array<double, 4> fv, fvp;

    assert(f.size() > 0);

    c.resize(f.size());

    for ( auto g = 0; g < f.size(); ++g )
    {
        c[g].resize(extents[a[1].nx][a[0].nx][4][4]);
        for ( auto i = 0; i < a[1].nx; ++i )
        {
            for ( auto j = 0; j < a[0].nx; ++j )
            {
                for ( auto k = 0; k < 4; ++k )
                {
                    for ( auto l = 0; l < 4; ++l )
                    {
                        fv[l] = get_value_clamped(g,  i + k - 1, j + l - 1);
                    }
                    fvp = mmul(fv);
                    for ( auto l = 0; l < 4; ++l )
                    {
                        c[g][i][j][l][k] = fv[l];
                    }
                }
            }
        }
    }
}

inline std::array<double, 4>
bicubic_interpolator::mmul(const std::array<double, 4>& x)
{
    std::array<double, 4> rval;
    std::fill( rval.begin(), rval.end(), 0.0 );

    static double m[4][4] = { {0, 2, 0, 0}, {-1, 0, 1, 0}, {2, -5, 4, -1}, {-1, 3, -3, 1} };

    for ( auto i = 0; i < 4; ++i )
        for ( auto j = 0; j < 4; ++j )
            rval[i] += m [i][j] * x[i];

    return rval;
}

double
bicubic_interpolator::interp_val(size_t idx, double x, double y)
{

    auto const itx = std::lower_bound ( a[0].x.begin(), a[0].x.end(), x ) - 1;
    auto const ity = std::lower_bound ( a[1].x.begin(), a[1].x.end(), y ) - 1;

    auto ix = std::distance ( a[0].x.begin(), itx );
    auto iy = std::distance ( a[1].x.begin(), ity );

    CLAMP(ix, 0, a[0].nx - 1);
    CLAMP(iy, 0, a[1].nx - 1);

    double xp = ( x - a[0].x[ix] ) * a[0].idx;
    double yp = ( y - a[1].x[iy] ) * a[1].idx;


    std::array<double, 4> bv, bvp;
    for ( auto i = 0 ; i < 4; i++)
        bv[i] =0.5 *( c[idx][iy][ix][i][0]
                      + c[idx][iy][ix][i][1] * xp
                      + c[idx][iy][ix][i][2] * xp * xp
                      + c[idx][iy][ix][i][3] * xp * xp * xp);


    bvp = mmul(bv);

    double rval = 0.5 *( bvp[0]
                         + bvp[1] * yp
                         + bvp[2] * yp * yp
                         + bvp[3] * yp * yp * yp);

    return rval;
}


bicubic_interpolator::bicubic_interpolator ( const axes2D &input_axis, const std::vector<data2D> &input_data )
{
    set_data( input_axis, input_data );
    make_c();
}

bicubic_interpolator::bicubic_interpolator ( const std::string &ncdf_file )
{
    load_from_netcdf ( ncdf_file );
    make_c();
}


double
bicubic_interpolator::operator () (size_t idx, double x, double y)
{
    return interp_val(idx, x, y);
}

void
bicubic_interpolator::set_data ( const axes2D &input_axis, const std::vector<data2D> &input_data )
{
    a = input_axis;
    f.assign ( input_data.begin(), input_data.end() );
}

void
bicubic_interpolator::load_from_netcdf ( const std::string &nc_file )
{
    std::vector<std::string> dims = {"temperature", "saturation"};
    std::vector<std::string> cols = {"nucleation_rate", "critical_size"};

    netcdf_reader nr ( nc_file, dims, cols );

    set_data ( nr.get_dims(), nr.get_data() );
}
