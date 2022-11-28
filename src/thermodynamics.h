#pragma once
#include <cmath>
#include "data_structures.h"
#include "initial_cond.h"
#include <iostream>

namespace thermo
{

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

        Y[0] = (rho1/rho);
        return Y;
    }

    inline double r_mix(std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double r = 0;

        for(int i = 0; i < mass_frac.size(); i++)
        {
            r += mass_frac[i]*par[i].r;
        }
        
        return r;
    }

    inline double gamma_mix(std::vector<double> const& mass_frac, std::vector<parameters> const& par)
    {
        double gamma = 0;
        
        for(int i = 0; i < mass_frac.size(); i++)
        {
            gamma += mass_frac[i]*par[i].gamma;
        }
        return gamma;
    }

    inline double pressure(int vel_comp,int n_comp, std::vector<parameters> const& par,double* U)
    {
        double p = 0;

        for(int i = n_comp; i < n_comp + vel_comp; i++)
        {
            p += U[i]*U[i];
        }

        double gamma = gamma_mix(composition(n_comp,U),par);

        return (gamma-1)*(U[n_comp+vel_comp]-0.5*p/U[0]);
    }

    inline double temperature(int vel_comp,int n_comp, std::vector<parameters> const& par, double* U)
    {
        double p = pressure(vel_comp,n_comp,par,U);
        double r = r_mix(composition(n_comp,U),par);
        return p/U[0]/r;
    }

    inline double speed_of_sound(int vel_comp, int n_comp, std::vector<parameters> const& par, double*U)
    {
        double p = pressure(vel_comp,n_comp,par,U);
        double gamma = gamma_mix(composition(n_comp,U),par);
        //std::cout << composition(n_comp,U)[0] << " " << composition(n_comp,U)[1] << " " << composition(n_comp,U)[2] << "\n";
        //std::cout << p << "\n";
        return sqrt(gamma*p/U[0]);
    }

    inline double r_mix(initial_conditions const& IC, std::vector<parameters> const& par)
    {
        double r = 0;

        for(int i = 0; i < IC.composition_mass_frac.size(); i++)
        {
            r += IC.composition_mass_frac[i]*par[i].r;
        }
        return r;
    }

    inline double gamma_mix(initial_conditions const& IC, std::vector<parameters> const& par) // ????
    {
        double gamma = 0;

        for(int i = 0; i < IC.composition_mass_frac.size(); i++)
        {
            gamma += IC.composition_mass_frac[i]*par[i].gamma;
        }
        return gamma;
    }

    inline double mach_number(int vel_comp,int n_comp, std::vector<parameters> const& par, double* U)
    {   
        double U2 = 0;

        for(int i = 0; i < vel_comp; i++)
        {
            U2 += pow(U[n_comp+i]/U[0],2); 
        }

        double c = speed_of_sound(vel_comp,n_comp,par,U);

        return sqrt(U2)/c;
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

