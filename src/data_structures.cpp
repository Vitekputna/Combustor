#include "data_structures.h"
#include "thermodynamics.h"
#include "mesh.h"
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <iostream>
#include <cstring>

typedef unsigned int uint;

parameters::parameters(){}

array::array(){}

array::array(int N, int k) : k{k}
{
    allocate(N,k);
}

void array::allocate(int N, int m_k)
{
    k = m_k;
    this->N = N;
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

array& array::operator=(array& other)
{
    k = other.k;
    N = other.N;
    std::memcpy(arr,other.arr,N);

    return *this;
}

variables::variables(int N, int N_walls, int dim, int vel_comp, int n_comp, int N_res) : N{N}, dim{dim}, N_walls{N_walls}, N_res{N_res}, vel_comp{vel_comp}, n_comp{n_comp}
{
    W.allocate(N*dim,dim);
    wall_flux.allocate(N_walls*dim,dim);
    diff_flux.allocate(N_walls*dim,dim);
    grad.allocate(2*N*dim,2*dim);
    alfa.allocate(N*dim,dim);
    Source.allocate(N*dim,dim);
    wall_grad.allocate(N*dim,dim);  

    T_grad = (double*)(calloc(N_walls,sizeof(double)));
    p = (double*)(malloc(N*sizeof(double)));
    T = (double*)(malloc(N*sizeof(double)));
    M = (double*)(malloc(N*sizeof(double)));
    res = (double*)(calloc((int)(N_res),sizeof(double)));
}

variables::variables(int N, int N_walls, int dim, int vel_comp, int n_comp, int N_res, std::vector<double>& U) : variables(N, N_walls, dim, vel_comp, n_comp, N_res)
{
    for(uint n = 0; n < N; n++)
    {  
        for(uint i = 0; i < dim; i++)
        {
            W(n,i) = U[i];
        }
    }
}

variables::variables(mesh const msh, config const cfg, std::vector<std::vector<double>> U) : variables(msh.N, msh.N_walls, cfg.dim, cfg.vel_comp, cfg.n_comp, cfg.max_iter/cfg.n_r)
{
    for(auto const& group : msh.physical_surface)
    {
        // std::cout << group.group_value << "\n";
        for(auto const& idx : group.member_idx)
        {
            // std::cout << idx << "\n";
            for(uint i = 0; i < dim; i++)
            {
                // std::cout << i << "\n";
                // std::cout << U.size() << "\n";
                W(idx,i) = U[group.group_value][i];
            }
        }
    }
}

variables::~variables()
{
    free(T_grad);
    free(p);
    free(T);
    free(M);
    free(res);
}

void variables::pressure(std::vector<parameters> const& par)
{
    for(uint i = 0; i < N; i++)
    {
        p[i] = thermo::pressure(vel_comp,n_comp,par,W(i));
    }
}

void variables::temperature(std::vector<parameters> const& par)
{
    for(uint i = 0; i < N; i++)
    {
        T[i] = thermo::temperature(vel_comp,n_comp,par,W(i));
    }
}

void variables::mach_number(std::vector<parameters> const& par)
{
    for(uint i = 0; i < N; i++)
    {
        M[i] = thermo::mach_number(vel_comp,n_comp,par,W(i));
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

double variables::compute_T_vertex_average(vertex const& node, mesh const& msh)
{
    double T_avg;

    for(auto const& cell_idx : node.common_cell_idx)        
    {
        T_avg += (1/(double)node.common_cell_idx.size())*T[cell_idx];
    }

    return T_avg;
}