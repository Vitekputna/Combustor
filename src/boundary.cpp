#include <iostream>
#include "boundary.h"
#include "mesh.h"

boundary::boundary(mesh const& msh) : msh{msh} {}

void boundary::apply(variables& var)
{
    for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    {
        BC_funcs[msh.ghost_cell_val[i]](msh.ghost_cell_idx[i],var,
                                        msh,bc_val + msh.ghost_cell_val[i]*4);
    }
}

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

void boundary::subsonic_inlet(int idx, variables& var, mesh const& msh, double* P)
{
    std::cout << "sub inlet : " << idx << "\n";
}

void boundary::subsonic_outlet(int idx, variables& var, mesh const& msh, double* P)
{
    std::cout << "sub outlet : " << idx << "\n";
}