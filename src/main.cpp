#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"
#include "vtk_export.h"
#include "solver.h"
#include "initial_cond.h"
#include "file_input.h"

#include "startup.h"

int N_threads;

int main(int argc, char** argv)
{
    N_threads = std::stoi(argv[2]);
    std::cout << message;
    parameters par;
    config cfg;
    initial_conditions IC;
    solver sol;

    mesh msh("mesh/" + std::string(argv[1]));
    
    boundary bdr(msh,par,cfg);
    
    read_config_files(msh, bdr, cfg, par, IC,sol);
    
    variables var(msh.N,msh.N_walls,cfg.dim,cfg.vel_comp,cfg.max_iter/cfg.n_r,IC.U);
    
    sol.solve(var,msh,bdr,par,cfg);

    export_vtk(var,msh,"exp.vtk");
    export_res(var, "res.txt");
}   