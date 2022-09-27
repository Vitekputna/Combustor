#pragma once
#include "data_structures.h"
#include "mesh.h"

void compute_wall_flux(variables& var, mesh const& msh, parameters const& par, void(*flux)(int,int,int,
                                                                                           variables&,
                                                                                           parameters const&,
                                                                                           face const&))
{
    int n,o;

    #pragma omp parallel for private(o,n) shared(msh,var,par)
        for(int w = 0; w < msh.N_walls;w++)
        {
            n = msh.walls[w].neigbour_cell_index;
            o = msh.walls[w].owner_cell_index;

            flux(w,n,o,var,par,msh.walls[w]);
        }
}

void compute_cell_res(variables& var, mesh const& msh, config const& cfg)
{
    int f;
    #pragma omp parallel for private(f) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        for(uint k = 0; k < var.dim; k++)
        {
            f = 0;
            for(auto const& wall : msh.cells[c].cell_walls)
            {
                var.W(c,k) -= cfg.dt/msh.cells[c].V*var.wall_flux(wall,k)*msh.cells[c].owner_idx[f];
                f++;
            }
        }
    }
}

void compute_cell_gradient(variables& var, mesh const& msh)
{
    
}