#include <iostream>
#include "boundary.h"
#include "mesh.h"

boundary::boundary(mesh const& msh) : msh{msh} {}

void boundary::apply(variables& var, double* bc_val)
{


    for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    {
        BC_funcs[msh.ghost_cell_val[i]](msh.ghost_cell_idx[i],var,
                                        msh,bc_val + msh.ghost_cell_val[i]*4);
    }
}

void boundary::wall(int idx, variables& var,mesh const& msh, double* P) //wall symmetry
{
    int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

    double nx = msh.walls[msh.cells[idx].cell_walls[0]].n[0];
    double ny = msh.walls[msh.cells[idx].cell_walls[0]].n[1];

    var.rho[idx] = var.rho[cell_idx];
    var.e[idx] = var.e[cell_idx];

    var.rhou[idx] = var.rhou[cell_idx] - 2*(var.rhou[cell_idx]*nx + var.rhov[cell_idx]*ny)*nx;
    var.rhov[idx] = var.rhov[cell_idx] - 2*(var.rhou[cell_idx]*nx + var.rhov[cell_idx]*ny)*ny;
}

void boundary::supersonic_inlet(int idx, variables& var, mesh const& msh, double* P) // set p0,T0, u magintude, u direction
{
    var.rho[idx] = P[0];
    var.rhou[idx] = P[1];
    var.rhov[idx] = P[2];
    var.e[idx] = P[3];
}

void boundary::supersonic_outlet(int idx, variables& var, mesh const& msh, double* P) // copy all
{
    int cell_idx = msh.walls[msh.cells[idx].cell_walls[0]].owner_cell_index;

    var.rho[idx] = var.rho[cell_idx];
    var.e[idx] = var.e[cell_idx];
    var.rhou[idx] = var.rhou[cell_idx];
    var.rhov[idx] = var.rhov[cell_idx];
}

void boundary::subsonic_inlet(int idx, variables& var, mesh const& msh, double* P)
{
    std::cout << "sub inlet : " << idx << "\n";
}

void boundary::subsonic_outlet(int idx, variables& var, mesh const& msh, double* P)
{
    std::cout << "sub outlet : " << idx << "\n";
}