#include "data_structures.h"
#include "boundary.h"
#include "mesh.h"
#include "numerical_flux.h"
#include "time_step.h"

#include <limits>
#include <algorithm>

typedef unsigned int uint;

void solve(variables& var, mesh& msh, boundary& bdr, parameters& par, config& cfg,double* bc_val)
{
    std::cout << "////////////////////////////////////////////////////\n";
    std::cout << "Computation running...\n\n";

    int n,o;
    int f;
    int t = 1;
    int r = 0;
    double delta;
    double res = cfg.max_res*2;

    //for(uint t = 1; t < cfg.iter; t++) 
    do
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
                f = 0;
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    var.W(c,k) -= cfg.dt/msh.cells[c].V*var.wall_flux(wall,k)*msh.cells[c].owner_idx[f];
                    f++;
                }
            }
        }

        if(!(t % cfg.n_b))
        {
            bdr.apply(var,bc_val);
        }
        
        if(!(t % cfg.n_t))
        {
            cfg.dt = cfg.CFL*time_step(msh,par,var);
        }

        if(!(t % cfg.n_r))
        {
            res = 0;
            for(uint c = 0; c < msh.N_cells; c++)
            {
                f = 0,delta = 0;
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    delta += var.wall_flux(wall,cfg.res_idx)*msh.cells[c].owner_idx[f];
                    f++;
                }
                
                res = std::max(res,abs(delta));
            }

            std::cout << "                                       \r"; 
            std::cout << "Time iteration: " <<  t << "\t"
                      << "Residual: " <<  res << "\t\r" << std::flush;

            var.res[r] = res;
            r++;
        }
        t++;
        
    } while((res > cfg.max_res || t < cfg.min_iter) && t < cfg.max_iter);
    
    std::cout << "\n";

    var.pressure(par); 
    var.temperature(par);
    var.mach_number(par);
}