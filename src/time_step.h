#pragma once
#include "data_structures.h"
#include "mesh.h"
#include "thermodynamics.h"
#include <limits.h>

double time_step(mesh const& msh,parameters const& par,config const& cfg,variables& var)
{
    double u = 0,c;
    double dt =  std::numeric_limits<double>::max();
    for(int i = 0; i < var.N; i++)
    {
        for(int j = 0; j < var.vel_comp; j++)
        {
            u += var.W(i,var.n_comp+j)*var.W(i,var.n_comp+j);
        }

        u = sqrt(u);

        c = sqrt(par.gamma*thermo::pressure(var.vel_comp,var.n_comp,par,var.W(i))/var.W(i,0));

        dt = std::min(dt,msh.min_V/(2*(c+u)));
    }

    return dt;
}
