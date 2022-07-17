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

    double U[4] = {0.1,1,2,5.5};

    variables var(msh.N,msh.N_walls,4,U);
    parameters par;
    boundary bdr(msh);

    double P[20] = {0,0,0,0,
                   1,2,3,4,
                   0,0,0,0,
                   0,0,0,0,
                   0,0,0,0};

    bdr.apply(var,P);
    
    for(uint i = 0; i < var.N; i++)
    {
        for(uint k = 0; k < var.dim; k++)
        {
            std::cout << var.W(i,k) << " ";
        }
        std::cout << "\n";
    }

    solve(var,msh,bdr,par,P);


    export_vtk(var,msh,"exp.vtk");

}