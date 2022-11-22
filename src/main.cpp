#include <iostream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <vector>
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
    // std::cout << message;
    std::vector<parameters> par; // specie parameters
    config cfg;
    solver sol;

    mesh msh("mesh/" + std::string(argv[1]));

    boundary bdr(msh,par,cfg);
    
    std::vector<std::vector<double>> IC_vec;
    IC_vec = read_config_files(msh, bdr, cfg, par, sol);

    variables var(msh, cfg, IC_vec);

    //sol.solve(var,msh,bdr,par,cfg);

    bdr.apply(var);

    var.pressure(par);
    var.temperature(par);
    var.mach_number(par);
    export_vtk(var,msh,"exp.vtk");
    export_res(var, "res.txt");
}
