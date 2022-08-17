#include <iostream>
#include <cmath>
#include "thermodynamics.h"
#include "boundary.h"
#include "mesh.h"
#include "nonlinear_solver.h"

typedef unsigned int uint;

boundary::boundary(mesh const& msh, parameters const& par, config& cfg) : msh{msh}, par{par}, cfg{cfg} {}

void boundary::apply(variables& var)
{
    // for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    // {
    //     (this->*BC_funcs[msh.ghost_cell_val[i]])(msh.ghost_cell_idx[i],var,
    //                                              msh,bc_val + msh.ghost_cell_val[i]*4);
    // }

    for(auto const& group : msh.boundary_groups)
    {
        (this->*BC_funcs[group.group_value])(group.member_idx,var,
                                                msh,bc_val + group.group_value*4);   
    }
}



void boundary::apply(variables& var, double* bc_val)
{
    // for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    // {
    //     (this->*BC_funcs[msh.ghost_cell_val[i]])(msh.ghost_cell_idx[i],var,
    //                                              msh,bc_val + msh.ghost_cell_val[i]*4);
    // }

    // for(auto const& group : msh.boundary_groups)
    // {
    //     for(uint i = 0; i < group.member_idx.size(); i++)
    //     {
    //         (this->*BC_funcs[group.group_value])(group.member_idx[i],var,
    //                                              msh,bc_val + group.group_value*4);
    //     }
    // }

    for(auto const& group : msh.boundary_groups)
    {   
        (this->*BC_funcs[group.group_value])(group.member_idx,var,
                                                msh,bc_val + group.group_value*4);
    }
}

void boundary::wall(std::vector<uint> group_idx, variables& var,mesh const& msh, double* P) //wall symmetry
{
    int cell_idx;
    double nx, ny;

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        nx = msh.walls[msh.cells[idx].cell_walls[0]].n[0];
        ny = msh.walls[msh.cells[idx].cell_walls[0]].n[1];

        var.W(idx,0) = var.W(cell_idx,0);
        var.W(idx,3) = var.W(cell_idx,3);

        var.W(idx,1) = var.W(cell_idx,1) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*nx;
        var.W(idx,2) = var.W(cell_idx,2) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*ny;
    }
}

void boundary::supersonic_inlet(std::vector<uint> group_idx, variables& var, mesh const& msh, double* P) // set p0,T0, u magintude, u direction
{
    for(auto const& idx : group_idx)
    {
        for(uint i = 0; i < var.dim; i++)
        {
            var.W(idx,i) = P[i];
        }
    }
}

void boundary::supersonic_outlet(std::vector<uint> group_idx, variables& var, mesh const& msh, double* P) // copy all
{
    int cell_idx;
    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        for(uint i = 0; i < var.dim; i++)
        {
            var.W(idx,i) = var.W(cell_idx,i);
        }
    }
}

inline double M_iter_func(boundary const& B,double M, double* P)
{
    return (P[1]/(B.par.gamma-1) + B.par.gamma*P[1]/2*M*M) * pow(1+(B.par.gamma-1)/2*M*M,B.par.gamma/(1-B.par.gamma))-P[0];
}

void boundary::subsonic_inlet(std::vector<uint> group_idx, variables& var, mesh const& msh, double* P)
{
    int cell_idx;
    double e = 0 ,Min,c;

    int N = group_idx.size();

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;
        e += var.W(cell_idx,3);
    }
    e = e/N;

    //compute inlet mach number
    double params[5] = {0,1,cfg.bisec_iter*1.0,e,P[0]};
    Min = bisection_method(M_iter_func, *this, params, 2);

    double rho = thermo::isoentropic_density(par,par.r*P[1]/P[0],Min); //density
    c = sqrt(par.gamma*par.r*thermo::isoentropic_temperature(par,P[1],Min));
    double ru = rho*c*Min*cos(P[2]);
    double rv = rho*c*Min*sin(P[2]);

    for(auto const& idx : group_idx)
    {
        var.W(idx,0) = rho;
        var.W(idx,1) = ru;
        var.W(idx,2) = rv;
        var.W(idx,3) = e;
    }
}

void boundary::subsonic_outlet(std::vector<uint> group_idx, variables& var, mesh const& msh, double* P)
{
    int cell_idx;

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        var.W(idx,0) = var.W(cell_idx,0);
        var.W(idx,1) = var.W(cell_idx,1);
        var.W(idx,2) = var.W(cell_idx,2);

        var.W(idx,3) = P[0]/(par.gamma-1) + 0.5*(var.W(idx,1)+var.W(idx,2))/var.W(idx,0);
    }

    
}

// void boundary::wall(int idx, variables& var,mesh const& msh, double* P) //wall symmetry
// {
//     int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

//     double nx = msh.walls[msh.cells[idx].cell_walls[0]].n[0];
//     double ny = msh.walls[msh.cells[idx].cell_walls[0]].n[1];

//     var.W(idx,0) = var.W(cell_idx,0);
//     var.W(idx,3) = var.W(cell_idx,3);

//     var.W(idx,1) = var.W(cell_idx,1) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*nx;
//     var.W(idx,2) = var.W(cell_idx,2) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*ny;
// }

// void boundary::supersonic_inlet(int idx, variables& var, mesh const& msh, double* P) // set p0,T0, u magintude, u direction
// {
//     for(uint i = 0; i < var.dim; i++)
//     {
//         var.W(idx,i) = P[i];
//     }
// }

// void boundary::supersonic_outlet(int idx, variables& var, mesh const& msh, double* P) // copy all
// {
//     int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

//     for(uint i = 0; i < var.dim; i++)
//     {
//         var.W(idx,i) = var.W(cell_idx,i);
//     }
// }

// inline double M_iter_func(boundary const& B,double M, double* P)
// {
//     return (P[1]/(B.par.gamma-1) + B.par.gamma*P[1]/2*M*M) * pow(1+(B.par.gamma-1)/2*M*M,B.par.gamma/(1-B.par.gamma))-P[0];
// }

// void boundary::subsonic_inlet(int idx, variables& var, mesh const& msh, double* P)
// {
//     int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

//     double e = var.W(cell_idx,3);

//     //compute inlet mach number
//     double params[5] = {0,1,cfg.bisec_iter*1.0,e,P[0]};
//     double Min = bisection_method(M_iter_func, *this, params, 2);

//     var.W(idx,0) = thermo::isoentropic_density(par,par.r*P[1]/P[0],Min); //density

//     double c = sqrt(par.gamma*par.r*thermo::isoentropic_temperature(par,P[1],Min));

//     var.W(idx,1) = var.W(idx,0)*c*Min*cos(P[2]);
//     var.W(idx,2) = var.W(idx,0)*c*Min*sin(P[2]);
//     var.W(idx,3) = e;
// }

// void boundary::subsonic_outlet(int idx, variables& var, mesh const& msh, double* P)
// {
//     int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

//     var.W(idx,0) = var.W(cell_idx,0);
//     var.W(idx,1) = var.W(cell_idx,1);
//     var.W(idx,2) = var.W(cell_idx,2);

//     var.W(idx,3) = P[0]/(par.gamma-1) + 0.5*(var.W(idx,1)+var.W(idx,2))/var.W(idx,0);
// }