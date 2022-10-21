#pragma once
#include <cmath>
#include "data_structures.h"

namespace thermo
{
    inline double pressure(int dim, parameters const& par,double* U)
    {
        return (par.gamma-1)*(U[dim-1]-0.5*(U[1]*U[1] + U[2]*U[2] + (dim > 2)*U[3]*U[3])/U[0]);
    }

    inline double temperature(int dim, parameters const& par,double* U)
    {
        double p = pressure(dim,par,U);
        return p/U[0]/par.r;
    }

    inline double mach_number(int dim, parameters const& par, double* U)
    {   
        double p = pressure(dim,par,U);
        double k = par.gamma;
        double t = (-1+sqrt(1-4*(k-1)/k*(1/(k-1)-U[dim-1]/p)))/(k-1);

        return sqrt(t);
    }

    inline double mach_number_stagnate(int dim, parameters const& par, double* U)
    {
        double T = temperature(dim,par,U);
        double M = sqrt((U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/U[0]/U[0])/sqrt(par.gamma*par.r*T);

        return sqrt((-1+sqrt(1+2*(par.gamma-1)*M*M))/(par.gamma-1));
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

