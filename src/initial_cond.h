#include "data_structures.h"
#include "cmath"

void supersonic_inlet(double* U,double p_0,double T_0, double Min, double alfa,parameters par)
{
    U[4] = pow(1+(par.gamma-1)/2*Min*Min,1/(1-par.gamma))*p_0/par.r/T_0;
    double p = pow(1+(par.gamma-1)/2*Min*Min,par.gamma/(1-par.gamma))*p_0;
    double c = sqrt(par.gamma*p/U[4]);
    U[5] = c*Min*cos(alfa)*U[4];
    U[6] = c*Min*sin(alfa)*U[4];
    U[7] = p/(par.gamma-1) + 0.5*U[4]*(U[5]*U[5] + U[6]*U[6]);
}

void no_move_flow(double* U,double p_0,double T_0,parameters par)
{
    U[0] = p_0/par.r/T_0;
    U[1] = 0;
    U[2] = 0;
    U[3] = p_0/(par.gamma-1);
}