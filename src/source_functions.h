#pragma once
#include "data_structures.h"
#include "thermodynamics.h"
#include "mesh.h"

inline void no_source_cartesian(variables& var, mesh const& msh, config const& cfg, parameters const& par){}

inline void heat_source_cartesian(variables& var, mesh const& msh, config const& cfg, parameters const& par)
{
    int f;
    #pragma omp parallel for private(f) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        var.W(c,var.dim-1) += cfg.dt*1e5;
    }
}

inline void axisymetric_source(variables& var, mesh const& msh, config const& cfg, parameters const& par)
{
    int f;
    #pragma omp parallel for private(f) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        var.W(c,2) += cfg.dt*var.p[c]/msh.cells[c].y;
        var.W(c,2) += cfg.dt*var.W(c,0)*var.W(c,3)*var.W(c,3)/msh.cells[c].y;
        var.W(c,3) += -cfg.dt*var.W(c,0)*var.W(c,2)*var.W(c,3)/msh.cells[c].y;
    }
}
