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

int N_threads;

int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));
    N_threads = std::stoi(argv[2]);

    msh.export_mesh();

    parameters par;

    config cfg;
    cfg.n_t = 100;
    cfg.n_r = 1000;
    cfg.max_res = 1;
    cfg.CFL = 1;
    cfg.bisec_iter = 15;
    cfg.n_r = 5000;
    cfg.max_iter = 3e6;

    boundary bdr(msh,par,cfg);

    initial_conditions IC;
    IC.alfa = 0;
    IC.Min = 1.5;
    IC.M_start = 1.5;
    IC.p_0 = 1e5;
    IC.p_start = 1e5;
    IC.p_stat =7.5e4;
    IC.T_0 = 300;
    IC.par = par;

    no_move_flow(IC);

    variables var(msh.N,msh.N_walls,4,cfg.max_iter/cfg.n_r,IC.U);
    
    supersonic_inlet(IC);
    supersonic_outlet(IC);
    //subsonic_inlet(P,1e5,300,0,par);
    //subsonic_outlet(P,7.5e4,par);

    bdr.apply(var,IC.B);

    solve(var,msh,bdr,par,cfg,IC.B);
    export_vtk(var,msh,"exp.vtk");
    export_res(var, "res.txt");



}   