#include <iostream>
#include <string>
#include "mesh.h"
#include "data_structures.h"
#include "boundary.h"

int main(int argc, char** argv)
{
    mesh msh("mesh/" + std::string(argv[1]));
    
    variables var(msh.N,4);
    parameters par;
    boundary bdr(msh);
    bdr.apply(var);

}