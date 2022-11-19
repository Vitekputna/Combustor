#include "initial_func.h"
#include "thermodynamics.h"
#include <iostream>

std::vector<double> move_flow(int vel_comp, int n_comp, initial_conditions& IC, std::vector<parameters> const& par)
{
    //hustoty chemických složek
    double r = thermo::r_mix(IC,par);
    double gamma = thermo::gamma_mix(IC,par);

    IC.U.push_back(IC.p_start/r/IC.T_start); // total density

    for(int n = 1; n < n_comp; n++)
    {
        IC.U.push_back(IC.composition_mass_frac[n]*IC.p_start/par[n].r/IC.T_start);
    }

    double p = IC.p_start;
    double p_o = pow(1+(gamma-1)/2*IC.M_start*IC.M_start,gamma/(gamma-1))*IC.p_start;
    double T_o = (1+(gamma-1)/2*IC.M_start*IC.M_start)*IC.T_start;
    double c = sqrt(gamma*r*T_o);

    switch (vel_comp)
    {
    case 2:
        IC.U.push_back(c*IC.M_start*cos(IC.alfa_start)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.alfa_start)*IC.U[0]);

        IC.U.push_back(p/(gamma-1) + 0.5*(IC.U[n_comp]*IC.U[n_comp] + IC.U[n_comp +1]*IC.U[n_comp +1])/IC.U[0]);

        return IC.U;
        break;

    case 3:
        IC.U.push_back(c*IC.M_start*cos(IC.alfa_start)*cos(IC.beta)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.alfa_start)*cos(IC.beta)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.beta)*IC.U[0]);

        IC.U.push_back(p/(gamma-1) + 0.5*(IC.U[n_comp]*IC.U[n_comp] + IC.U[n_comp + 1]*IC.U[n_comp + 1] + IC.U[n_comp + 2]*IC.U[n_comp + 2])/IC.U[0]);

        return IC.U;
        break;

    default:
        std::cout << "Error creating initial conditions...\n";
        exit(1);
    }
}

void supersonic_inlet(int vel_comp, int n_comp, std::vector<parameters> const& par, boundary_group& bdr)
{
    //hustoty chemických složek
    double r = thermo::r_mix(bdr.composition,bdr.composition_mass_frac,par);
    double gamma = thermo::gamma_mix(bdr.composition,bdr.composition_mass_frac,par);

    bdr.bc_val.push_back(bdr.p_0/r/bdr.T_0);

    for(int i = 1; i < n_comp; i++)
    {
        bdr.bc_val.push_back(bdr.composition_mass_frac[i]*bdr.p_0/par[i].r/bdr.T_0);
    }

    double p = bdr.p_0;
    double T_o = (1+(gamma-1)/2*bdr.Min*bdr.Min)*bdr.T_0;
    double c = sqrt(gamma*r*T_o);

    switch (vel_comp)
    {
    case 2:
        bdr.bc_val.push_back(c*bdr.Min*cos(bdr.alfa)*bdr.bc_val[0]);
        bdr.bc_val.push_back(c*bdr.Min*sin(bdr.alfa)*bdr.bc_val[0]);
        break;

    case 3: // pridat uhel rotace proudu
        bdr.bc_val.push_back(c*bdr.Min*cos(bdr.alfa)*cos(bdr.beta)*bdr.bc_val[0]);
        bdr.bc_val.push_back(c*bdr.Min*sin(bdr.alfa)*cos(bdr.beta)*bdr.bc_val[0]);
        bdr.bc_val.push_back(c*bdr.Min*sin(bdr.beta)*bdr.bc_val[0]);
        break;

    default:
        std::cout << "Error creating initial conditions...\n";
        exit(1);
    }

    bdr.bc_val.push_back(p/(gamma-1) + 0.5*(bdr.bc_val[1]*bdr.bc_val[1] + bdr.bc_val[2]*bdr.bc_val[2])/bdr.bc_val[0]);
}

void supersonic_outlet(std::vector<parameters> const& par, boundary_group& bdr)
{
    //nothing to do here
}

void subsonic_inlet(std::vector<parameters> const& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(bdr.p_0);
    bdr.bc_val.push_back(bdr.T_0);
    bdr.bc_val.push_back(bdr.alfa);
    bdr.bc_val.push_back(bdr.beta);
}

void subsonic_outlet(std::vector<parameters> const& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(bdr.p_stat);
}