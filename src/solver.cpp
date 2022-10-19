#include <iostream>
#include "solver.h"
#include "data_structures.h"
#include "boundary.h"
#include "mesh.h"
#include "time_step.h"
#include "fvm.h"
#include "vtk_export.h"
#include "math.h"

#include <limits>
#include <algorithm>

#include "omp.h"

extern int N_threads;

typedef unsigned int uint;

solver::solver(){}

void solver::solve(variables& var, mesh& msh, boundary& bdr, parameters& par, config& cfg)
{
    std::cout << "////////////////////////////////////////////////////\n";
    std::cout << "Computation running...\n\n";

    omp_set_num_threads(N_threads);
    std::cout << "Running on: " << omp_get_max_threads() << " threads\n";

    int f;
    int t = 1;
    int r = 0;
    double delta; 
    double time = 0;
    double last_time = 0;
    double res = cfg.max_res*2;

    do 
    {  
        var.pressure(par);
        //compute_cell_gradient(var,msh);
        //grad_limiting(var,msh);
        compute_wall_flux(var,msh,par,flux_func);
        compute_cell_res(var,msh,cfg,par,source_func);

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
            std::cout << "Time : " <<  time << " s\t"
                      << "Residual: " <<  res << "\t\r" << std::flush;

            var.res[r] = res;
            r++;
        }

        if(time - last_time >= cfg.export_interval)        
        {
            last_time = time;
            var.pressure(par); 
            var.temperature(par);
            var.mach_number(par);

            //std::cout << "out/" + msh.name.substr(5,msh.name.length()-5) + "_" + std::to_string(t) + ".vtk" << "\n";
            
            export_vtk(var,msh,"timesteps/" + msh.name.substr(5,msh.name.length()-5) + "_" + std::to_string(time) + ".vtk");
        }

        t++;
        time += cfg.dt;
        
   } while((res > cfg.max_res || t < cfg.min_iter) && t < cfg.max_iter && time < cfg.max_time);
    
    std::cout << "\n";

    var.temperature(par);
    var.mach_number(par);
}