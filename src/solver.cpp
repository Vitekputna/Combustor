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

#include <omp.h>

extern int N_threads;

typedef unsigned int uint;

solver::solver(){}

double inline solver::runtime_checks(int t, variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg)
{
    double res = cfg.max_res*2;

    if(!(t % cfg.n_b))
    {
        bdr.apply(var);
    }
    
    if(!(t % cfg.n_t))
    {
        cfg.dt = cfg.CFL*time_step(msh,par,cfg,var);
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
        var.mach_number(par);
        var.pressure(par);
        var.temperature(par);
        export_vtk(var,msh,"timesteps/" + msh.name.substr(5,msh.name.length()-5) + "_" + std::to_string(time) + ".vtk");
    }

    return res;
}

void solver::explicit_euler(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg)
{
    double res = cfg.max_res*2;

    do
    {  
        compute_wall_flux(var,msh,par,flux_func);
        apply_cell_res(var,msh,cfg,par,source_func);

        res = runtime_checks(t,var,msh,bdr,par,cfg);

        t++;
        time += cfg.dt;

    } while((res > cfg.max_res || t < cfg.min_iter) && t < cfg.max_iter && time < cfg.max_time);

    std::cout << "Computation completed, final residue: " << res << "\n"; 
}

void solver::Adams_Bashforth(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg)
{
    double res = cfg.max_res*2;

    array old_res;
    old_res.allocate(var.N*var.dim,var.dim);

    array new_res;
    new_res.allocate(var.N*var.dim,var.dim);

    compute_wall_flux(var,msh,par,flux_func);
    compute_cell_res(old_res,var,msh,cfg,par,source_func);
    apply_cell_res(var,msh,cfg,par,source_func);

    do
    {  
        compute_wall_flux(var,msh,par,flux_func);
        compute_cell_res(new_res,var,msh,cfg,par,source_func);
        apply_cell_res(old_res,new_res,-0.5,1.5,var,msh,cfg,par);

        old_res = new_res;

        res = runtime_checks(t,var,msh,bdr,par,cfg);

        t++;
        time += cfg.dt;

    } while((res > cfg.max_res || t < cfg.min_iter) && t < cfg.max_iter && time < cfg.max_time);

    std::cout << "Computation completed, final residue: " << res << "\n"; 
}

void solver::solve(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg)
{
    std::cout << "Number of species: " << var.n_comp << ", Number of velocity components: " << var.vel_comp << "\n";
    std::cout << "////////////////////////////////////////////////////\n";
    std::cout << "Computation running...\n\n";

    bdr.apply(var);

    Adams_Bashforth(var,msh,bdr,par,cfg);

    var.pressure(par);
    var.temperature(par);
    var.mach_number(par);
}