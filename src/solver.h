#include "data_structures.h"
#include "boundary.h"
#include "mesh.h"
#include "numerical_flux.h"
#include "time_step.h"

void solve(variables& var, mesh& msh, boundary& bdr, parameters& par, config& cfg,double* bc_val)
{
    int n,o;

    for(uint t = 0; t < cfg.iter; t++)
    {
        for(uint w = 0; w < msh.N_walls;w++)
        {
            n = msh.walls[w].neigbour_cell_index;
            o = msh.walls[w].owner_cell_index;

            HLL_flux(w,n,o,var,par,msh.walls[w]);
        }

        for(uint c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                int f = 0;
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    var.W(c,k) -= cfg.dt/msh.cells[c].V*var.wall_flux(wall,k)*msh.cells[c].owner_idx[f];
                    f++;
                }
            }
        }

        bdr.apply(var,bc_val);

        if(!(t % cfg.n_t))
        {
            cfg.dt = time_step(var);
        }

    }
}