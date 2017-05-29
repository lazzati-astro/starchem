#include <cmath>
#include <plog/Log.h>

#include "reaction.h"

double
reaction::rate ( double tgas )
{
    double k = 0.0;
    switch ( type )
    {
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 14:
        k = alpha * pow ( tgas / 300.0, beta ) * exp ( -gamma / tgas );
        break;
    default:
        LOGE << "UNKNOWN TYPE (" << type << ")";
    }
    return k;
}


