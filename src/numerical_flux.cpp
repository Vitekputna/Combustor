#include "numerical_flux.h"
#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "math.h"
#include <iostream>

typedef unsigned int uint;

void HLL_flux(std::vector<uint> idx ,int vel_comp, int n_comp, double* w, double* n, double* o, std::vector<parameters> const& par, face const& f)
{
    int dim = n_comp+vel_comp+1;

    double So,Sn;
    double uo,vo;
    double po,pn;
    double un,vn;
    double co,cn; 

    double phi[dim], phi_o[dim], phi_n[dim];

    po = thermo::pressure(vel_comp,n_comp,par,o);
    pn = thermo::pressure(vel_comp,n_comp,par,n);

    double gamma_o = thermo::gamma_mix(thermo::composition(n_comp,o),par);
    double gamma_n = thermo::gamma_mix(thermo::composition(n_comp,n),par);

    co = sqrt(gamma_o*po/ o[0]);
    cn = sqrt(gamma_n*pn/ n[0]);

    uo = o[n_comp]/o[0];
    un = n[n_comp]/n[0];

    vo = o[n_comp+1]/o[0];
    vn = n[n_comp+1]/n[0];

    Sn = std::max(un*f.n[0]+vn*f.n[1] + cn,
                  uo*f.n[0]+vo*f.n[1] + co);

    So = std::min(un*f.n[0]+vn*f.n[1] - cn,
                  uo*f.n[0]+vo*f.n[1] - co);    

    double Xn[dim];
    double Xo[dim];

    for(int n = 0; n < n_comp; n++)
    {
        Xn[n] = 0;
        Xo[n] = 0;
    }

    double vel_X[3] = {f.n[0],f.n[1],0}; 

    for(int v = 0; v < vel_comp; v++)
    {
        Xn[v+n_comp] = vel_X[v];
        Xo[v+n_comp] = vel_X[v];
    }

    Xn[dim-1] = un*f.n[0]+vn*f.n[1];
    Xo[dim-1] = uo*f.n[0]+vo*f.n[1];

    // for(uint k = 0; k < dim; k++)
    for(auto const& k : idx)
    {
        phi_n[k] = (un*f.n[0]+vn*f.n[1])*n[k]+pn*Xn[k];
        phi_o[k] = (uo*f.n[0]+vo*f.n[1])*o[k]+po*Xo[k];
        phi[k] = (Sn*phi_o[k]-So*phi_n[k]+So*Sn*(n[k]-o[k]))/(Sn-So);

        if(Sn <= 0)
        {
            w[k] = phi_n[k]*f.S;
        }
        else if(So >= 0)
        {
            w[k] = phi_o[k]*f.S;
        }
        else
        {
            w[k] = phi[k]*f.S;
        }
    }
}

void Upwind_flux(std::vector<uint> idx ,int vel_comp, int n_comp, double* w, double* n, double* o, std::vector<parameters> const& par, face const& f)
{
    int dim = n_comp+vel_comp+1;

    double uo,vo;
    double un,vn;
    double uf,vf;

    double phi_o[dim], phi_n[dim];

    uo = o[n_comp]/o[0];
    un = n[n_comp]/n[0];

    vo = o[n_comp+1]/o[0];
    vn = n[n_comp+1]/n[0];

    uf = 0.5*(uo+un);
    vf = 0.5*(vo+vn);

    for(auto const& k : idx)
    {
        phi_n[k] = (un*f.n[0]+vn*f.n[1])*n[k];
        phi_o[k] = (uo*f.n[0]+vo*f.n[1])*o[k];

        if(uf*f.n[0]+vf*f.n[1] < 0)
        {
            w[k] = phi_n[k]*f.S;
        }
        else if (uf*f.n[0]+vf*f.n[1] > 0)
        {
            w[k] = phi_o[k]*f.S;
        }
        else
        {
            w[k] = 0;
        }
    }
}