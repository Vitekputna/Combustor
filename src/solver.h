#pragma once
#include "data_structures.h"
#include "source_functions.h"
#include "numerical_flux.h"
#include "mesh.h"
#include "boundary.h"


class solver
{
public:

    int f;
    int t = 1;
    int r = 0;
    double delta; 
    double time = 0;
    double last_time = 0;

    void(*flux_func)(std::vector<uint>,int,int,double*,double*,double*,std::vector<parameters> const&,face const&) = HLL_flux;
    void(*source_func)(array&,variables&,mesh const&,config const&,std::vector<parameters> const&) = ternary_chem;
    solver();

    void solve(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg);

    double inline runtime_checks(int t, variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg);
    void explicit_euler(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg);
    void Adams_Bashforth(variables& var, mesh& msh, boundary& bdr, std::vector<parameters>& par, config& cfg);
};