#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "math.h"

typedef unsigned int uint;

void HLL_flux(int w, int n, int o, variables & var,parameters const& par, face const& f)
{
    double So,Sn;
    double uo,vo;
    double po,pn;
    double un,vn;
    double co,cn;

    double phi[4], phi_o[4], phi_n[4];

    po = thermo::pressure(par,&var.W(o,0));
    pn = thermo::pressure(par,&var.W(n,0));

    co = sqrt(par.gamma*po/var.W(o,0));
    cn = sqrt(par.gamma*pn/var.W(n,0));

    uo = var.W(o,1)/var.W(o,0);
    un = var.W(n,1)/var.W(n,0);

    vo = var.W(o,2)/var.W(o,0);
    vn = var.W(n,2)/var.W(n,0);

    Sn = std::max(un*f.n[0]+vn*f.n[1] + cn,
                  uo*f.n[0]+vo*f.n[1] + co);

    So = std::min(un*f.n[0]+vn*f.n[1] - cn,
                  uo*f.n[0]+vo*f.n[1] - co);    

    double Xn[4] = {0,f.n[0],f.n[1],un*f.n[0]+vn*f.n[1]};    
    double Xo[4] = {0,f.n[0],f.n[1],uo*f.n[0]+vo*f.n[1]};    

    for(uint k = 0; k < var.dim; k++)
    {
        phi_n[k] = (un*f.n[0]+vn*f.n[1])*var.W(n,k)+pn*Xn[k];
        phi_o[k] = (uo*f.n[0]+vo*f.n[1])*var.W(o,k)+po*Xo[k];
        phi[k] = (Sn*phi_o[k]-So*phi_n[k]+So*Sn*(var.W(n,k)-var.W(o,k)))/(Sn-So);

        if(Sn <= 0)
        {
            var.wall_flux(w,k) = phi_n[k]*f.S;
        }
        else if(So >= 0)
        {
            var.wall_flux(w,k) = phi_o[k]*f.S;
        }
        else
        {
            var.wall_flux(w,k) = phi[k]*f.S;
        }
    } 
}