#include <iostream>
#include "boundary.h"

boundary::boundary(mesh const& msh) : msh{msh} {}

void boundary::apply(variables& var)
{

    for(uint i = 0; i < msh.ghost_cell_idx.size(); i++)
    {
        BC_funcs[msh.ghost_cell_val[i]](msh.ghost_cell_idx[i],var);
    }

}

void boundary::wall(int idx, variables& var)
{
    std::cout << "wall : " << idx << " " << var.rho[idx] << "\n";
}

void boundary::supersonic_inlet(int idx, variables& var)
{
    std::cout << "sup inlet : " << idx << "\n";
}

void boundary::supersonic_outlet(int idx, variables& var)
{
    std::cout << "sup outlet : " << idx << "\n";
}

void boundary::subsonic_inlet(int idx, variables& var)
{
    std::cout << "sub inlet : " << idx << "\n";
}

void boundary::subsonic_outlet(int idx, variables& var)
{
    std::cout << "sub outlet : " << idx << "\n";
}