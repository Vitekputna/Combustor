#include <iostream>
#include <cmath>
#include "thermodynamics.h"
#include "boundary.h"
#include "mesh.h"

boundary::boundary(mesh const& msh, parameters const& par) : msh{msh}, par{par} {}

void boundary::apply(variables& var)
{
    for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    {
        (this->*BC_funcs[msh.ghost_cell_val[i]])(msh.ghost_cell_idx[i],var,
                                                 msh,bc_val + msh.ghost_cell_val[i]*4);

        // BC_funcs[msh.ghost_cell_val[i]](msh.ghost_cell_idx[i],var,
        //                                 msh,bc_val + msh.ghost_cell_val[i]*4);
    }
}

void boundary::apply(variables& var, double* bc_val)
{
    for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    {
        (this->*BC_funcs[msh.ghost_cell_val[i]])(msh.ghost_cell_idx[i],var,
                                                 msh,bc_val + msh.ghost_cell_val[i]*4);
    }
}

void boundary::wall(int idx, variables& var,mesh const& msh, double* P) //wall symmetry
{
    int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

    double nx = msh.walls[msh.cells[idx].cell_walls[0]].n[0];
    double ny = msh.walls[msh.cells[idx].cell_walls[0]].n[1];

    var.W(idx,0) = var.W(cell_idx,0);
    var.W(idx,3) = var.W(cell_idx,3);

    var.W(idx,1) = var.W(cell_idx,1) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*nx;
    var.W(idx,2) = var.W(cell_idx,2) - 2*(var.W(cell_idx,1)*nx + var.W(cell_idx,2)*ny)*ny;
}

void boundary::supersonic_inlet(int idx, variables& var, mesh const& msh, double* P) // set p0,T0, u magintude, u direction
{
    for(uint i = 0; i < var.dim; i++)
    {
        var.W(idx,i) = P[i];
    }
}

void boundary::supersonic_outlet(int idx, variables& var, mesh const& msh, double* P) // copy all
{
    int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

    for(uint i = 0; i < var.dim; i++)
    {
        var.W(idx,i) = var.W(cell_idx,i);
    }
}

inline double boundary::M_iter_func(double M, double e, double* P)
{
    return (P[0]/(par.gamma-1) + par.gamma*P[0]/2*M*M) * pow(1+(par.gamma-1)/2*M*M,par.gamma/(1-par.gamma))-e;
}

void boundary::subsonic_inlet(int idx, variables& var, mesh const& msh, double* P)
{

    std::cout << "called\n";

    int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

    double e = var.W(cell_idx,3);

    //compute inlet mach number
    double M1 = 0;
    double dM = 0.01;
    double M2 = M1 + dM;

    double Min;

    while (M2 <= 1+dM)
    {
        if(M_iter_func(M1,e,P)*M_iter_func(M2,e,P) <= 0)
        {
            std::cout << M1 << " " << M2 << "\n";

            Min = 0.5*(M1+M2);

            break;
        }

        M1 += dM;
        M2 += dM;
    }

    var.W(idx,0) = thermo::isoentropic_density(par,par.r*P[1]/P[0],Min); //density

    double c = sqrt(par.gamma*par.r*thermo::isoentropic_temperature(par,P[1],Min));

    var.W(idx,1) = var.W(idx,0)*c*Min*cos(P[2]);
    var.W(idx,2) = var.W(idx,0)*c*Min*sin(P[2]);
    var.W(idx,3) = e;
}

void boundary::subsonic_outlet(int idx, variables& var, mesh const& msh, double* P)
{
    
}