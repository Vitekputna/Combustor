#pragma once
#include "data_structures.h"
#include "boundary.h"
#include "cmath"

struct initial_conditions
{
    double U[4]; //initial cond
    double B[20]; //boundary value cond
    double p_0, T_0, Min, alfa, p_stat, p_start, M_start, T_start, alfa_start;
    parameters par;
};

void no_move_flow(initial_conditions& IC)
{
    IC.U[0] = IC.p_start/IC.par.r/IC.T_start;
    IC.U[1] = 0;
    IC.U[2] = 0;
    IC.U[3] = IC.p_start/(IC.par.gamma-1);
}

void move_flow(initial_conditions& IC)
{
    IC.U[0] = pow(1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start,1/(1-IC.par.gamma))*IC.p_start/IC.par.r/IC.T_start;
    double p = pow(1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start,IC.par.gamma/(1-IC.par.gamma))*IC.p_start;
    double c = sqrt(IC.par.gamma*p/IC.U[0]);
    IC.U[1] = c*IC.M_start*cos(IC.alfa_start)*IC.U[0];
    IC.U[2] = c*IC.M_start*sin(IC.alfa_start)*IC.U[0];
    IC.U[3] = p/(IC.par.gamma-1) + 0.5*(IC.U[1]*IC.U[1] + IC.U[2]*IC.U[2])/IC.U[0];
}

void supersonic_inlet(parameters& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(pow(1+(par.gamma-1)/2*bdr.Min*bdr.Min,1/(1-par.gamma))*bdr.p_0/par.r/bdr.T_0);
    double p = pow(1+(par.gamma-1)/2*bdr.Min*bdr.Min,par.gamma/(1-par.gamma))*bdr.p_0;
    double c = sqrt(par.gamma*par.r*bdr.T_0);
    bdr.bc_val.push_back(c*bdr.Min*cos(bdr.alfa)*bdr.bc_val[0]);
    bdr.bc_val.push_back(c*bdr.Min*sin(bdr.alfa)*bdr.bc_val[0]);
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
}

void subsonic_outlet(parameters& par, boundary_group& bdr)
{
    bdr.bc_val.push_back(bdr.p_stat);
}