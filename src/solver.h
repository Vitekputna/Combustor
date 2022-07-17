#include <iostream>
#include "data_structures.h"
#include "boundary.h"
#include "mesh.h"


void solve(variables& var, mesh& msh, boundary& bdr, parameters& par,double* bc_val)
{
    int n,o;

    for(uint t = 0; t < 2; t++)
    {
        for(uint w = 0; w < msh.N_walls;w++)
        {
            n = msh.walls[w].neigbour_cell_index;
            o = msh.walls[w].owner_cell_index;

            for(uint k = 0; k < var.dim; k++)
            {
                var.wall_flux(w,k);
            }
        }

        for(uint c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    var.W(c,k) -= 1e-1/msh.cells[c].V*var.wall_flux(wall,k);
                }
            }
        }
    }

    bdr.apply(var,bc_val);
}