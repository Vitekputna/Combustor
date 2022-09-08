#pragma once
#include "data_structures.h"
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
    std::cout << IC.M_start << "\n";

    IC.U[0] = pow(1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start,1/(1-IC.par.gamma))*IC.p_start/IC.par.r/IC.T_start;
    double p = pow(1+(IC.par.gamma-1)/2*IC.M_start*IC.M_start,IC.par.gamma/(1-IC.par.gamma))*IC.p_start;
    double c = sqrt(IC.par.gamma*p/IC.U[0]);
    IC.U[1] = c*IC.M_start*cos(IC.alfa_start)*IC.U[0];
    IC.U[2] = c*IC.M_start*sin(IC.alfa_start)*IC.U[0];
    IC.U[3] = p/(IC.par.gamma-1) + 0.5*(IC.U[1]*IC.U[1] + IC.U[2]*IC.U[2])/IC.U[0];

    std::cout << IC.U[0] << " " << IC.U[1] << " " << IC.U[2] << " " << IC.U[3] << "\n";
    std::cout << c << " " << p << "\n";
}

void supersonic_inlet(initial_conditions& IC)
{
    IC.B[4] = pow(1+(IC.par.gamma-1)/2*IC.Min*IC.Min,1/(1-IC.par.gamma))*IC.p_0/IC.par.r/IC.T_0;
    double p = pow(1+(IC.par.gamma-1)/2*IC.Min*IC.Min,IC.par.gamma/(1-IC.par.gamma))*IC.p_0;
    double c = sqrt(IC.par.gamma*IC.par.r*IC.T_0);
    IC.B[5] = c*IC.Min*cos(IC.alfa)*IC.B[4];
    IC.B[6] = c*IC.Min*sin(IC.alfa)*IC.B[4];
    IC.B[7] = p/(IC.par.gamma-1) + 0.5*(IC.B[5]*IC.B[5] + IC.B[6]*IC.B[6])/IC.B[4];   
}

void supersonic_outlet(initial_conditions& IC)
{
    //nothing to do here
}

void subsonic_inlet(initial_conditions& IC)
{
    IC.B[12] = IC.p_0;
    IC.B[13] = IC.T_0;
    IC.B[14] = IC.alfa;
    IC.B[15] = 0;
}

void subsonic_outlet(initial_conditions& IC)
{
   IC.B[16] = IC.p_stat;
}