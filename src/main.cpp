#include <iostream>
#include <string>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"
#include "vtk_export.h"
#include "solver.h"
#include "initial_cond.h"

#include "thermodynamics.h"

int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));
    //msh.export_mesh();

    parameters par;
    boundary bdr(msh);
    config cfg;
    cfg.dt = 1e-6;
    cfg.iter = 1e3;
    cfg.n_t = 100;

    //initial
    double U[4];
    no_move_flow(U,1e5,300,par);
    variables var(msh.N,msh.N_walls,4,U);

    //boundary
    double P[20];
    supersonic_inlet(P,1e5,300,1.5,0,par);

    solve(var,msh,bdr,par,cfg,P);
    export_vtk(var,msh,"exp.vtk");
}