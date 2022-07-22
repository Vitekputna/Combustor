#include <iostream>
#include <string>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"
#include "vtk_export.h"
#include "solver.h"


int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));

    double U[4] = {1.225,0,0,250000};

    variables var(msh.N,msh.N_walls,4,U);
    parameters par;
    boundary bdr(msh);
    config cfg;
    cfg.dt = 1e-4;

    double P[20] = {0,0,0,0,
                   1.225,40,0,400000,
                   0,0,0,0,
                   0,0,0,0,
                   0,0,0,0};

    solve(var,msh,bdr,par,cfg,P);
    export_vtk(var,msh,"exp.vtk");

}