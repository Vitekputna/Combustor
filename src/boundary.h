#pragma once
#include "mesh.h"
#include "data_structures.h"
#include <vector>

class boundary
{
    typedef void (boundary::*func)(int,variables&,mesh const&,double*);

    

    public:
    boundary(mesh const& msh,parameters const& par);

    mesh const& msh;
    parameters const& par;

    double bc_val[20];

    void apply(variables& var); 
    void apply(variables& var, double* bc_val); 

    void wall(int idx, variables& var, mesh const& msh, double* P); //0
    void supersonic_inlet(int idx, variables& var, mesh const& msh, double* P); //1
    void supersonic_outlet(int idx, variables& var, mesh const& msh, double* P); //2
    void subsonic_inlet(int idx, variables& var, mesh const& msh, double* P); //3
    void subsonic_outlet(int idx, variables& var, mesh const& msh, double* P); //4

    // std::vector<func> BC_funcs = {&boundary::wall};

    std::vector<func> BC_funcs = {&boundary::wall,&boundary::supersonic_inlet
                                 ,&boundary::supersonic_outlet,&boundary::subsonic_inlet
                                 ,&boundary::subsonic_outlet};

    friend inline double M_iter_func(boundary const& B, double M, double* P);                                 
};

