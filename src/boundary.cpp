#include <iostream>
#include <cmath>
#include <algorithm>
#include "thermodynamics.h"
#include "boundary.h"
#include "mesh.h"
#include "nonlinear_solver.h"
#include "math.h"

typedef unsigned int uint;

boundary::boundary(mesh const& msh, parameters const& par, config& cfg) : msh{msh}, par{par}, cfg{cfg} {}


void boundary::apply(variables& var)
{
    int g = 0;
    for(auto const& group : boundary_groups)
    {   
        (this->*BC_funcs[group.bc_func_idx])(group.member_idx,var,msh,group.bc_val);
        g++;
    }
}

void boundary::wall(std::vector<uint> const& group_idx, variables& var,mesh const& msh, std::vector<double> const& P) //wall symmetry
{
    int cell_idx;
    double nx, ny;

    std::vector<double> n(3);

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        n[0] = msh.walls[msh.cells[idx].cell_walls[0]].n[0];
        n[1] = msh.walls[msh.cells[idx].cell_walls[0]].n[1];
        n[2] = 0;

        var.W(idx,0) = var.W(cell_idx,0);
        var.W(idx,var.dim-1) = var.W(cell_idx,var.dim-1);

        for(int i = 0; i < var.vel_comp; i++)
        {
            var.W(idx,i+1) = var.W(cell_idx,i+1) - 2*(var.W(cell_idx,1)*n[0] + var.W(cell_idx,2)*n[1])*n[i];
        }

        //var.W(idx,1) = var.W(cell_idx,1) - 2*(var.W(cell_idx,1)*n[0] + var.W(cell_idx,2)*n[1])*n[0];
        //var.W(idx,2) = var.W(cell_idx,2) - 2*(var.W(cell_idx,1)*n[0] + var.W(cell_idx,2)*n[1])*n[1];
    }
}

void boundary::supersonic_inlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P) // set p0,T0, u magintude, u direction
{
    for(auto const& idx : group_idx)
    {
        for(uint i = 0; i < var.dim; i++)
        {
            var.W(idx,i) = P[i];
        }
    }
}

void boundary::supersonic_outlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P) // copy all
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

void boundary::subsonic_inlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P)
{
    int cell_idx;
    double e = 0 ,Min,c;
    double p;
    double ru_b, rv_b, rw_b;
    double u_b, v_b, w_b;
    double rho_b;

    int N = group_idx.size();

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;
        e += var.W(cell_idx,var.dim-1);
        p += var.p[cell_idx];

        ru_b = var.W(cell_idx,1);
        rv_b = var.W(cell_idx,2);
        rw_b = var.W(cell_idx,3);

        rho_b += var.W(cell_idx,0);
    }
    e = e/N;
    p = p/N;
    ru_b = ru_b/N;
    rv_b = rv_b/N;
    rw_b = rw_b/N;
    rho_b = rho_b/N;

    u_b = ru_b/rho_b;
    v_b = rv_b/rho_b;
    w_b = rw_b/rho_b;


    //compute inlet mach number
    // double params[5] = {0,1,cfg.bisec_iter*1.0,e,P[0]};
    // Min = bisection_method(M_iter_func, *this, params, 2);

    //Min = std::max(0.0,sqrt((-1+sqrt(1-4*(-1)/par.gamma*(1/(par.gamma-1)-e/p)))/(par.gamma-1)));

    //std::cout << e << " " << p << " " << Min <<"\n";

    double rho = P[0]/par.r/P[1]; //density
    //double T_o = (1+(par.gamma-1)/2*Min*Min)*P[1];
    //c = sqrt(par.gamma*par.r*P[1]);
    // double ru = rho*c*Min*cos(P[2])*cos(P[3]);
    // double rv = rho*c*Min*sin(P[2])*cos(P[3]);
    // double rw = rho*c*Min*sin(P[3]);


    double e_boundary = P[0]/(par.gamma-1) + 0.5*rho_b*(u_b*u_b + v_b*v_b + w_b*w_b);

    std::vector<double> V = {u_b,v_b,w_b};

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        var.W(idx,0) = rho;

        for(int i = 0; i < var.vel_comp; i++)
        {
            //var.W(idx,i+1) = V[i];
            var.W(idx,i+1) = var.W(cell_idx,i+1);
        }

        var.W(idx,var.dim-1) = e_boundary;
    }
}

void boundary::subsonic_outlet(std::vector<uint> const& group_idx, variables& var, mesh const& msh, std::vector<double> const& P)
{
    int cell_idx;
    //double M_max;

    // for(auto const& idx : group_idx)
    // {   
    //     cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;
    //     M_max = std::max(0.0,std::abs(thermo::mach_number(par,var.W(cell_idx))));
    // }

    for(auto const& idx : group_idx)
    {
        cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

        var.W(idx,0) = var.W(cell_idx,0);
        // var.W(idx,1) = var.W(cell_idx,1);
        // var.W(idx,2) = var.W(cell_idx,2);
    
        for(int i = 0; i < var.vel_comp; i++)
        {
            var.W(idx,i+1) = var.W(cell_idx,i+1);
        }

        var.W(idx,var.dim-1) = P[0]/(par.gamma-1) + 0.5*(var.W(idx,1)*var.W(idx,1)+var.W(idx,2)*var.W(idx,2))/var.W(idx,0);

        //std::cout << var.W(idx,var.dim-1) << "\n";

        // if(M_max >= 1)
        // {
        //     var.W(idx,3) = var.W(cell_idx,3);
        // }
        // else
        // {
        //     var.W(idx,3) = P[0]/(par.gamma-1) + 0.5*(var.W(idx,1)+var.W(idx,2))/var.W(idx,0);
        // }
    }
}

boundary_group::boundary_group(){}

boundary_group::boundary_group(int i) : bc_func_idx{i} {}

boundary_group::boundary_group(int i, const std::vector<unsigned int> idxs) : bc_func_idx{i}, member_idx{idxs} {}