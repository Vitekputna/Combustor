#include "initial_func.h"
#include <iostream>

void no_move_flow(int dim, initial_conditions& IC)
{
    int vel_comp = dim-2;

    IC.U.push_back(IC.p_start/IC.par.r/IC.T_start);
    
    for(int i = 0; i < vel_comp; i++) {IC.U.push_back(0);}

    IC.U.push_back(IC.p_start/(IC.par.gamma-1));
}

std::vector<double> move_flow(int vel_comp, int n_comp, initial_conditions& IC)
{
    //hustoty chemických složek
    for(int n = 0; n < n_comp; n++)
    {
        IC.U.push_back(IC.p_start/IC.par.r/IC.T_start);
    }

    double p = IC.p_start;
    double p_o = pow(1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start,IC.par.gamma/(IC.par.gamma-1))*IC.p_start;
    double T_o = (1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start)*IC.T_start;
    double c = sqrt(IC.par.gamma*IC.par.r*T_o);

    switch (vel_comp)
    {
    case 2:
        IC.U.push_back(c*IC.M_start*cos(IC.alfa_start)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.alfa_start)*IC.U[0]);

        IC.U.push_back(p/(IC.par.gamma-1) + 0.5*(IC.U[n_comp]*IC.U[n_comp] + IC.U[n_comp +1]*IC.U[n_comp +1])/IC.U[0]);

        return IC.U;
        break;

    case 3:
        IC.U.push_back(c*IC.M_start*cos(IC.alfa_start)*cos(IC.beta)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.alfa_start)*cos(IC.beta)*IC.U[0]);
        IC.U.push_back(c*IC.M_start*sin(IC.beta)*IC.U[0]);

        IC.U.push_back(p/(IC.par.gamma-1) + 0.5*(IC.U[n_comp]*IC.U[n_comp] + IC.U[n_comp + 1]*IC.U[n_comp + 1] + IC.U[n_comp + 2]*IC.U[n_comp + 2])/IC.U[0]);

        return IC.U;
        break;

    default:
        std::cout << "Error creating initial conditions...\n";
        exit(1);
    }
}

void supersonic_inlet(int vel_comp, int n_comp, parameters& par, boundary_group& bdr)
{
    for(int i = 0; i < n_comp; i++)
    {
        bdr.bc_val.push_back(bdr.p_0/par.r/bdr.T_0);
    }

    double p = bdr.p_0;
    double T_o = (1+(par.gamma-1)/2*bdr.Min*bdr.Min)*bdr.T_0;
    double c = sqrt(par.gamma*par.r*T_o);

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

    bdr.bc_val.push_back(p/(par.gamma-1) + 0.5*(bdr.bc_val[1]*bdr.bc_val[1] + bdr.bc_val[2]*bdr.bc_val[2])/bdr.bc_val[0]);
}

void supersonic_outlet(parameters& par, boundary_group& bdr)
{
    //nothing to do here
}

void subsonic_inlet(parameters& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(bdr.p_0);
    bdr.bc_val.push_back(bdr.T_0);
    bdr.bc_val.push_back(bdr.alfa);
    bdr.bc_val.push_back(bdr.beta);
}

void subsonic_outlet(parameters& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(bdr.p_stat);
}