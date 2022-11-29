#pragma once
#include "data_structures.h"
#include "source_functions.h"
#include "numerical_flux.h"
#include "mesh.h"
#include "boundary.h"


class solver
{
public:
    void(*flux_func)(int,int,double*,double*,double*,std::vector<parameters> const&,face const&) = HLL_flux;
    void(*source_func)(variables&,mesh const&,config const&,std::vector<parameters> const&) = ternary_chem;
    solver();
    void solve(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg);
};