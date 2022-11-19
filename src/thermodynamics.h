#pragma once
#include <cmath>
#include "data_structures.h"
#include "initial_cond.h"
#include <iostream>

namespace thermo
{
    inline double pressure(int vel_comp,int n_comp, parameters const& par,double* U)
    {
        double p = 0;

        for(int i = n_comp; i < n_comp + vel_comp; i++)
        {
            p += U[i]*U[i];
        }

        return (par.gamma-1)*(U[n_comp+vel_comp]-0.5*p/U[0]);
    }

    inline double temperature(int vel_comp,int n_comp, parameters const& par, double* U)
    {
        double p = pressure(vel_comp,n_comp,par,U);
        return p/U[0]/par.r;
    }

    inline double speed_of_sound(int vel_comp, int n_comp, parameters const& par, double*U)
    {
        double p = pressure(vel_comp,n_comp,par,U);
        return sqrt(par.gamma*p/U[0]);
    }

    inline double r_mix(initial_conditions const& IC, std::vector<parameters> const& par)
    {
        double r = 0;

        int idx = 0;
        for(auto const& i : IC.composition)
        {
            r += IC.composition_mass_frac[idx]*par[i].r;
            idx++;
        }
        return r;
    }

    inline double r_mix(std::vector<int> const& composition, std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double r = 0;

        int idx = 0;
        for(auto const& i : composition)
        {
            r += mass_frac[idx]*par[i].r;
            idx++;
        }
        return r;
    }

    inline double r_mix(std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double r = 0;
        int idx = 0;
        for(auto const& spec : par)
        {
            r += mass_frac[idx]*spec.r;
        }
        return r;
    }

    inline double gamma_mix(initial_conditions const& IC, std::vector<parameters> const& par) // ????
    {
        double gamma = 0;

        int idx = 0;
        for(auto const& i : IC.composition)
        {
            gamma += IC.composition_mass_frac[idx]*par[i].gamma;
            idx++;
        }
        return gamma;
    }
    
    inline double gamma_mix(std::vector<int> const& composition, std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double gamma = 0;

        int idx = 0;
        for(auto const& i : composition)
        {
            gamma += mass_frac[idx]*par[i].r;
            idx++;
        }
        return gamma;
    }

    inline double gamma_mix(std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double gamma = 0;
        int idx = 0;
        for(auto const& spec : par)
        {
            gamma += mass_frac[idx]*spec.gamma;
        }
        return gamma;
    }

    inline std::vector<double> composition(int n_comp, double* W)
    {
        double rho = W[0];
        double rho1 = rho;

        std::vector<double> Y(n_comp,0.0);

        for(int i = 1; i < n_comp; i++)
        {
            Y[i] = W[i]/rho;
            rho1 -= W[i];
        }

        Y[0] = (rho1/rho)*(rho1/rho > 0);

        return Y;
    }

    inline double mach_number(int vel_comp,int n_comp, parameters const& par, double* U)
    {   
        double p = pressure(vel_comp,n_comp,par,U);
        double k = par.gamma;
        double t = (-1+sqrt(1-4*(k-1)/k*(1/(k-1)-U[n_comp+vel_comp]/p)))/(k-1);

        return sqrt(t);
    }

    // inline double mach_number_stagnate(int dim, parameters const& par, double* U)
    // {
    //     double T = temperature(dim,par,U);
    //     double M = sqrt((U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/U[0]/U[0])/sqrt(par.gamma*par.r*T);

    //     return sqrt((-1+sqrt(1+2*(par.gamma-1)*M*M))/(par.gamma-1));
    // }

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

