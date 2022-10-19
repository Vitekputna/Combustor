#pragma once
#include "data_structures.h"
#include "source_functions.h"
#include "numerical_flux.h"
#include "mesh.h"
#include "boundary.h"


class solver
{
public:
    void(*flux_func)(int,double*,double*,double*,parameters const&,face const&) = HLL_flux_axi;
    void(*source_func)(variables&,mesh const&,config const&,parameters const&) = no_source_cartesian;

    solver();

    void solve(variables& var, mesh& msh, boundary& bdr, parameters& par, config& cfg);
};