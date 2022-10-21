#pragma once
#include "mesh.h"
#include "data_structures.h"
#include <vector>

typedef unsigned int uint;

struct boundary_group
{   
    boundary_group();
    boundary_group(int i);
    boundary_group(int i, const std::vector<unsigned int> idxs);
    std::vector<unsigned int> member_idx;
    int bc_func_idx;
    double p_0, T_0, Min, alfa, beta, p_stat;
    std::vector<double> bc_val;
};

class boundary
{
    typedef void (boundary::*func)(std::vector<uint> const&,variables&,mesh const&,std::vector<double> const&);

    public:
    boundary(mesh const& msh,parameters const& par, config& cfg);
    boundary();

    mesh const& msh;
    parameters const& par;
    config& cfg;

    double bc_val[20];

    void apply(variables& var); 

    void wall(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P); //0
    void supersonic_inlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P); //1
    void supersonic_outlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P); //2
    void subsonic_inlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P); //3
    void subsonic_outlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P); //4

    std::vector<func> BC_funcs = {&boundary::wall,
                                  &boundary::supersonic_inlet,
                                  &boundary::supersonic_outlet,
                                  &boundary::subsonic_inlet,
                                  &boundary::subsonic_outlet};

    std::vector<int> boundary_func_mask;

    std::vector<boundary_group> boundary_groups;

    friend inline double M_iter_func(boundary const& B, double M, double* P);                                 
};

