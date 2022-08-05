#pragma once
#include <cmath>
#include "data_structures.h"

namespace thermo
{
    inline double pressure(parameters const& par,double* U)
    {
        return (par.gamma-1)*(U[3]-0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);
    }

    inline double temperature(parameters const& par,double* U)
    {
        double p = pressure(par,U);
        return p/U[0]/par.r;
    }

    inline double isoentropic_pressure(parameters const& par, double p_0, double M)
    {
        return pow(1+(par.gamma-1)/2*M*M,par.gamma/(1-par.gamma))*p_0;
    }

    inline double isoentropic_density(parameters const& par, double rho_0, double M)
    {
        return pow(1+(par.gamma-1)/2*M*M,1/(1-par.gamma))*rho_0;
    }

    inline double isoentropic_temperature(parameters const& par, double T_0, double M)
    {
        return pow(1+(par.gamma-1)/2*M*M,-1)*T_0;
    }
}

