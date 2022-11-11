#pragma once
#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include <limits.h>

double time_step(mesh const& msh,parameters const& par,config const& cfg,variables& var)
{
    double u,c;
    double dt =  std::numeric_limits<double>::max();
    for(int i = 0; i < var.N; i++)
    {
        u = sqrt(var.W(i,1)*var.W(i,1) + var.W(i,2)*var.W(i,2));
        c = sqrt(par.gamma*thermo::pressure(var.vel_comp,var.n_comp,par,var.W(i))/var.W(i,0));

        dt = std::min(dt,msh.min_V/(2*(c+u)));
    }

    return cfg.CFL*dt;
}
