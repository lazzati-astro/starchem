#pragma once

#include <vector>
#include <array>
#include <boost/multi_array.hpp>

struct axis
{
    uint32_t nx;
    double   dx;
    double   idx;
    std::vector<double> x;

    void reset()
    {
        nx = 0;
        dx = 0;
        idx = 0;
        x.clear();
    }
};

typedef std::array<axis, 2>           axes2D;
typedef boost::multi_array<double, 2>        data2D;

