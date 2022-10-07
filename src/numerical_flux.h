#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "math.h"

typedef unsigned int uint;

void HLL_flux(int dim, double* w, double* n, double* o, parameters const& par, face const& f)
{
    double So,Sn;
    double uo,vo;
    double po,pn;
    double un,vn;
    double co,cn;

    double phi[4], phi_o[4], phi_n[4];

    po = thermo::pressure(par,o);
    pn = thermo::pressure(par,n);

    co = sqrt(par.gamma*po/ o[0]);
    cn = sqrt(par.gamma*pn/ n[0]);

    uo = o[1]/o[0];
    un = n[1]/n[0];

    vo = o[2]/o[0];
    vn = n[2]/n[0];

    Sn = std::max(un*f.n[0]+vn*f.n[1] + cn,
                  uo*f.n[0]+vo*f.n[1] + co);

    So = std::min(un*f.n[0]+vn*f.n[1] - cn,
                  uo*f.n[0]+vo*f.n[1] - co);    

    double Xn[4] = {0,f.n[0],f.n[1],un*f.n[0]+vn*f.n[1]};
    double Xo[4] = {0,f.n[0],f.n[1],uo*f.n[0]+vo*f.n[1]};

    for(uint k = 0; k < dim; k++)
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