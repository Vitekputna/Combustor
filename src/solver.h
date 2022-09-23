#include "data_structures.h"
#include "boundary.h"
#include "mesh.h"
#include "numerical_flux.h"
#include "time_step.h"
#include "fvm.h"

#include <limits>
#include <algorithm>

#include "omp.h"

extern int N_threads;

typedef unsigned int uint;

void solve(variables& var, mesh& msh, boundary& bdr, parameters& par, config& cfg,double* bc_val)
{
    std::cout << "////////////////////////////////////////////////////\n";
    std::cout << "Computation running...\n\n";

    omp_set_num_threads(N_threads);
    std::cout << "Running on: " << omp_get_max_threads() << " threads\n";

    int n,o;
    int f;
    int t = 1;
    int r = 0;
    double delta;
    double res = cfg.max_res*2;

    do 
    {
        compute_wall_flux(var,msh,par,HLL_flux);
        compute_cell_res(var,msh,cfg);
        

        if(!(t % cfg.n_b))
        {
            bdr.apply(var);
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

            std::cout << "                                                                         \r"; 
            std::cout << "Time iteration: " <<  t << "\t"
                      << "Residual: " <<  res << "\t\r" << std::flush;

            var.res[r] = res;
            r++;
        }

        if(!(t % cfg.n_exp))
        {
            var.pressure(par); 
            var.temperature(par);
            var.mach_number(par);

            //std::cout << "out/" + msh.name.substr(5,msh.name.length()-5) + "_" + std::to_string(t) + ".vtk" << "\n";
            
            export_vtk(var,msh,"timesteps/" + msh.name.substr(5,msh.name.length()-5) + "_" + std::to_string(t) + ".vtk");
        }

        t++;
        
   } while((res > cfg.max_res || t < cfg.min_iter) && t < cfg.max_iter);
    
    std::cout << "\n";

    var.pressure(par); 
    var.temperature(par);
    var.mach_number(par);
}