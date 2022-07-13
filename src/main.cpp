#include <iostream>
#include <string>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"

int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));

    double U[4] = {0.1,1,2,5.5};

    variables var(msh.N,4,U);
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
        std::cout << var.rho[i] << " " << var.rhou[i] << " " << var.rhov[i]
                  << " " << var.e[i] << "\n";
    }

}