#pragma once

#include <cmath>

namespace constants {
const double kBoltzmann = 1.38064852E-16;
const double pressure_0 = 1.0E6;

const double shapeFactor = pow(36.0 * M_PI, 1. / 3.);

double kBeta(const double T) {
    return kBoltzmann * T;
}
double equPres(const double A, const double B, const double T) {
    return exp( -A / T + B);
}

}
