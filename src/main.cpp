#include <iostream>
#include <string>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"
#include "vtk_export.h"
#include "solver.h"
#include "initial_cond.h"

int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));
    msh.export_mesh();

    parameters par;
    config cfg;
    boundary bdr(msh,par,cfg);
    cfg.n_t = 100;
    cfg.n_r = 1000;
    cfg.max_res = 1;
    cfg.CFL = 1;
    cfg.bisec_iter = 15;

    //initial
    double U[4];
    //no_move_flow(U,1e4,300,par);
    move_flow(U,1e4,300,0.2,0,par);
    variables var(msh.N,msh.N_walls,4,U);

    //boundary
    double P[20];
    supersonic_inlet(P,1e5,300,1.5,0,par);
    //subsonic_inlet(P,1e5,300,0,par);
    //subsonic_outlet(P,7.5e4,par);
    bdr.apply(var,P);

    solve(var,msh,bdr,par,cfg,P);
    export_vtk(var,msh,"exp.vtk");
    //export_res(var, "res.txt");



}   