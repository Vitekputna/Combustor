#include "data_structures.h"
#include "thermodynamics.h"
#include "mesh.h"
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <iostream>

typedef unsigned int uint;

array::array(){}

array::array(int N, int k) : k{k}
{
    allocate(N,k);
}

void array::allocate(int N, int m_k)
{
    k = m_k;
    arr = (double*)(calloc(N,sizeof(double)));
}

array::~array()
{
    free(arr);
}

double& array::operator()(int i, int j)
{
    return arr[i*k + j];
}

double* array::operator()(int i)
{
    return arr + i*k;
}

variables::variables(int N, int N_walls, int dim, int vel_comp, int N_res) : N{N}, dim{dim}, N_walls{N_walls}, N_res{N_res}, vel_comp{vel_comp}
{
    W.allocate(N*dim,dim);
    wall_flux.allocate(N_walls*dim,dim);
    grad.allocate(2*N*dim,2*dim);
    alfa.allocate(N*dim,dim);
    Source.allocate(N*dim,dim);

    p = (double*)(malloc(N*sizeof(double)));
    T = (double*)(malloc(N*sizeof(double)));
    M = (double*)(malloc(N*sizeof(double)));
    res = (double*)(calloc((int)(N_res),sizeof(double)));
}

variables::variables(int N, int N_walls, int dim, int vel_comp, int N_res, std::vector<double>& U) : variables(N, N_walls, dim, vel_comp, N_res)
{
    for(uint n = 0; n < N; n++)
    {  
        for(uint i = 0; i < dim; i++)
        {
            W(n,i) = U[i];
        }
    }
}

variables::variables(mesh const msh, config const cfg, std::vector<std::vector<double>> U) : variables(msh.N, msh.N_walls, cfg.dim, cfg.vel_comp, cfg.max_iter/cfg.n_r)
{
   for(auto const& group : msh.physical_surface)
   {
        for(auto const& idx : group.member_idx)
        {
            for(uint i = 0; i < dim; i++)
            {
                W(idx,i) = U[group.group_value][i];
            }
        }
   }
}

variables::~variables()
{
    free(p);
    free(T);
    free(M);
    free(res);
}

void variables::pressure(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        p[i] = thermo::pressure(dim,par,W(i));
    }
}

void variables::temperature(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        T[i] = thermo::temperature(dim,par,W(i));
    }
}

void variables::mach_number(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        M[i] = thermo::mach_number(dim,par,W(i));
    }
}

std::vector<double> variables::compute_vertex_average(vertex const& node, mesh const& msh)
{
    std::vector<double> W_avg(dim,0.0);

    for(int k = 0; k < dim; k++)
    {
        for(auto const& cell_idx : node.common_cell_idx)        
        {
            W_avg[k] += (1/(double)node.common_cell_idx.size())*W(cell_idx,k);
        }
    }

    return W_avg;
}