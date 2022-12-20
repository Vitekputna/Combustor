#pragma once
#include "data_structures.h"
#include "thermodynamics.h"
#include "mesh.h"

inline void no_source_cartesian(array& res, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par){}

inline void heat_source_cartesian(array& res, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    int f;
    #pragma omp parallel for private(f) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        var.W(c,var.dim-1) += cfg.dt*1e5;
    }
}

inline void axisymetric_source(array& res, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    int f;
    double p;

    #pragma omp parallel for private(f,p) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        p = thermo::pressure(cfg.vel_comp,cfg.n_comp,par,var.W(c));

        res(c,cfg.n_comp+1) += cfg.dt*p/msh.cells[c].y;
        res(c,cfg.n_comp+1) += cfg.dt*var.W(c,0)*var.W(c,cfg.n_comp+2)*var.W(c,cfg.n_comp+2)/msh.cells[c].y;
        res(c,cfg.n_comp+2) += -cfg.dt*var.W(c,0)*var.W(c,cfg.n_comp+1)*var.W(c,cfg.n_comp+2)/msh.cells[c].y;
    }
}

// inline void axisymetric_source(variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
// {
//     int f;
//     double p;

//     #pragma omp parallel for private(f,p) shared(var,cfg,msh)
//     for(int c = 0; c < msh.N_cells; c++)
//     {
//         p = thermo::pressure(var.vel_comp,var.n_comp,par,var.W(c));

//         var.W(c,var.n_comp+1) += cfg.dt*p/msh.cells[c].y;
//         var.W(c,var.n_comp+1) += cfg.dt*var.W(c,0)*var.W(c,var.n_comp+2)*var.W(c,var.n_comp+2)/msh.cells[c].y;
//         var.W(c,var.n_comp+2) += -cfg.dt*var.W(c,0)*var.W(c,var.n_comp+1)*var.W(c,var.n_comp+2)/msh.cells[c].y;
//     }
// }

inline void ternary_chem(array& res, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    double k = 1e5;
    std::vector<double> X;

    #pragma omp parallel for private(X) shared(var,cfg,msh,k)
    for(int c = 0; c < msh.N_cells; c++)
    {
        X = thermo::molar_composition(cfg.n_comp,var.W(c),par);

        res(c,0) += -cfg.dt*k*X[0]*X[1];
        res(c,1) += -cfg.dt*k*X[0]*X[1];
        res(c,2) += cfg.dt*k*X[0]*X[1];
    }   
}

// inline void ternary_chem(variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
// {
//     double k = 1e-3;
//     std::vector<double> X;

//     #pragma omp parallel for private(X) shared(var,cfg,msh,k)
//     for(int c = 0; c < msh.N_cells; c++)
//     {
//         X = thermo::molar_composition(var.n_comp,var.W(c),par);

//         var.W(c,0) += -cfg.dt*k*X[0]*X[1];
//         var.W(c,1) += -cfg.dt*k*X[0]*X[1];
//         var.W(c,2) += cfg.dt*k*X[0]*X[1];
//     }   
// }