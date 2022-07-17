#pragma once
#include "mesh.h"
#include "data_structures.h"
#include <vector>

typedef void (*func)(int,variables&,mesh const&,double*);

class boundary
{
    public:
    boundary(mesh const& msh);

    mesh const& msh;

    double bc_val[20];

    void apply(variables& var); 
    void apply(variables& var, double* bc_val); 

    static void wall(int idx, variables& var, mesh const& msh, double*); //0
    static void supersonic_inlet(int idx, variables& var, mesh const& msh, double*); //1
    static void supersonic_outlet(int idx, variables& var, mesh const& msh, double*); //2
    static void subsonic_inlet(int idx, variables& var, mesh const& msh, double*); //3
    static void subsonic_outlet(int idx, variables& var, mesh const& msh, double*); //4

    std::vector<func> BC_funcs = {&boundary::wall,&boundary::supersonic_inlet
                                 ,&boundary::supersonic_outlet,&boundary::subsonic_inlet
                                 ,&boundary::subsonic_outlet};
};