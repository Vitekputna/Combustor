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

int N_threads;

int main(int argc, char** argv)
{
    N_threads = std::stoi(argv[2]);

    parameters par;
    config cfg;
    initial_conditions IC;

    mesh msh("mesh/" + std::string(argv[1]));

    //msh.export_mesh();

    cfg.n_t = 100;
    cfg.max_res = 1;
    cfg.CFL = 1;
    cfg.bisec_iter = 15;
    cfg.n_r = 5000;
    cfg.max_iter = 1e7;
    cfg.n_b = 100;
    cfg.n_exp = 20000;

    boundary bdr(msh,par,cfg);

    read_config_files(bdr, cfg, par, IC);

    //IC.alfa = 0;
    //IC.Min = 1.5;
    //IC.M_start = 0.2;
    //IC.p_0 = 1e5;
    //IC.p_start = 7.5e4;
    //IC.p_stat =7.5e4;
    //IC.T_0 = 300;
    //IC.par = par;
    //IC.T_start = 300;
    //IC.alfa_start = 0;

    no_move_flow(IC);

    variables var(msh.N,msh.N_walls,4,cfg.max_iter/cfg.n_r,IC.U);
    
    supersonic_inlet(IC);
    supersonic_outlet(IC);
    subsonic_inlet(IC);
    subsonic_outlet(IC);

    bdr.apply(var,IC.B);

    solve(var,msh,bdr,par,cfg,IC.B);
    export_vtk(var,msh,"exp.vtk");
    export_res(var, "res.txt");



}   