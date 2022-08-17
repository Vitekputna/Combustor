#include "data_structures.h"
#include "cmath"

void no_move_flow(double* U,double p_0,double T_0,parameters par)
{
    U[0] = p_0/par.r/T_0;
    U[1] = 0;
    U[2] = 0;
    U[3] = p_0/(par.gamma-1);
}

void move_flow(double* U,double p_0,double T_0, double Min, double alfa,parameters par)
{
    U[0] = pow(1+(par.gamma-1)/2*Min*Min,1/(1-par.gamma))*p_0/par.r/T_0;
    double p = pow(1+(par.gamma-1)/2*Min*Min,par.gamma/(1-par.gamma))*p_0;
    double c = sqrt(par.gamma*p/U[0]);
    U[1] = c*Min*cos(alfa)*U[0];
    U[2] = c*Min*sin(alfa)*U[0];
    U[3] = p/(par.gamma-1) + 0.5*(U[1]*U[1] + U[2]*U[2])/U[0];
}

void supersonic_inlet(double* U,double p_0,double T_0, double Min, double alfa,parameters par)
{
    U[4] = pow(1+(par.gamma-1)/2*Min*Min,1/(1-par.gamma))*p_0/par.r/T_0;
    double p = pow(1+(par.gamma-1)/2*Min*Min,par.gamma/(1-par.gamma))*p_0;
    // double c = sqrt(par.gamma*p/U[4]);
    double c = sqrt(par.gamma*par.r*T_0);
    U[5] = c*Min*cos(alfa)*U[4];
    U[6] = c*Min*sin(alfa)*U[4];
    U[7] = p/(par.gamma-1) + 0.5*(U[5]*U[5] + U[6]*U[6])/U[4];
}

void subsonic_inlet(double* U,double p_0,double T_0, double alfa,parameters par)
{
    U[12] = p_0;
    U[13] = T_0;
    U[14] = alfa;
    U[15] = 0;
}

void subsonic_outlet(double* U,double p_stat,parameters par)
{
   U[16] = p_stat;
}