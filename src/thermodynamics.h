#pragma once
#include "data_structures.h"

inline double pressure(parameters const& par,double* U)
{
    return (par.gamma-1)*(U[3]-0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);
}

inline double temperature(parameters const& par,double* U)
{
    double p = pressure(par,U);
    return p/U[0]/par.r;
}