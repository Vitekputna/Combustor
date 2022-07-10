#pragma once
#include "mesh.h"
#include "data_structures.h"
#include <vector>

typedef void (*func)(int,variables&);

class boundary
{
    public:
    boundary(mesh const& msh);

    mesh const& msh;

    void apply(variables& var); 

    static void wall(int idx, variables& var); //0
    static void supersonic_inlet(int idx, variables& var); //1
    static void supersonic_outlet(int idx, variables& var); //2
    static void subsonic_inlet(int idx, variables& var); //3
    static void subsonic_outlet(int idx, variables& var); //4

    std::vector<func> BC_funcs = {&boundary::wall,&boundary::supersonic_inlet
                                 ,&boundary::supersonic_outlet,&boundary::subsonic_inlet
                                 ,&boundary::subsonic_outlet};

};


