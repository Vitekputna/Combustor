#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include "math.h"

void HLL_flux(int w, int n, int o, variables & var,parameters const& par, face const& f)
{
    double So,Sn;
    double uo,vo;
    double po,pn;
    double un,vn;
    double co,cn;

    double phi[4], phi_o[4], phi_n[4];

    //std::cout << "wall: " << w << " o: " << o << " n: " << n << "\n";

    po = pressure(par.gamma,&var.W(o,0));
    pn = pressure(par.gamma,&var.W(n,0));

    //std::cout << po << " " << pn << "\n";

    co = sqrt(par.gamma*po/var.W(o,0));
    cn = sqrt(par.gamma*pn/var.W(n,0));

    //std::cout << co << " " << cn << "\n";

    uo = var.W(o,1)/var.W(o,0);
    un = var.W(n,1)/var.W(n,0);

    vo = var.W(o,2)/var.W(o,0);
    vn = var.W(n,2)/var.W(n,0);

    //std::cout << uo << " " << un << "\n";
   // std::cout << vo << " " << vn << "\n";

    Sn = std::max(un*f.n[0]+vn*f.n[1] + cn,
                  uo*f.n[0]+vo*f.n[1] + co);

    So = std::min(un*f.n[0]+vn*f.n[1] - cn,
                  uo*f.n[0]+vo*f.n[1] - co);

    //std::cout << So << " " << Sn << "\n";
    

    double Xn[4] = {0,f.n[0],f.n[1],un*f.n[0]+vn*f.n[1]};    
    double Xo[4] = {0,f.n[0],f.n[1],uo*f.n[0]+vo*f.n[1]};    

    for(uint k = 0; k < var.dim; k++)
    {
        phi_n[k] = (un*f.n[0]+vn*f.n[1])*var.W(n,k)+pn*Xn[k];
        phi_o[k] = (uo*f.n[0]+vo*f.n[1])*var.W(o,k)+po*Xo[k];
        phi[k] = (Sn*phi_o[k]-So*phi_n[k]+So*Sn*(var.W(n,k)-var.W(o,k)))/(Sn-So);

        if(Sn <= 0)
        {
            //std::cout << "0: " << phi_n[k] << "\n";
            var.wall_flux(w,k) = phi_n[k]*f.S;
        }
        else if(So >= 0)
        {
            //std::cout<< "1: "  << phi_o[k] << "\n";
            var.wall_flux(w,k) = phi_o[k]*f.S;
        }
        else
        {
            //std::cout<< "2: "  << phi[k] << "\n";
            var.wall_flux(w,k) = phi[k]*f.S;
        }
    }
   // std::cout << phi_n[0] << " " << phi_o[0] << " " << phi[0] << "\n";
   // std::cout << phi_n[1] << " " << phi_o[1] << " " << phi[1] << "\n";
   // std::cout << phi_n[2] << " " << phi_o[2] << " " << phi[2] << "\n";
   // std::cout << phi_n[3] << " " << phi_o[3] << " " << phi[3] << "\n";
   // std::cout << "//////////////////////////////////\n";
}