#include "data_structures.h"
#include "thermodynamics.h"
#include "mesh.h"
#include <stdlib.h>
#include <limits>
#include <algorithm>

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

variables::variables(int N, int N_walls, int dim, int N_res) : N{N}, dim{dim}, N_walls{N_walls}, N_res{N_res}
{
    W.allocate(N*dim,dim);
    wall_flux.allocate(N_walls*dim,dim);
    grad.allocate(2*N*dim,2*dim);

    p = (double*)(malloc(N*sizeof(double)));
    T = (double*)(malloc(N*sizeof(double)));
    M = (double*)(malloc(N*sizeof(double)));
    res = (double*)(calloc((int)(N_res),sizeof(double)));
}

variables::variables(int N, int N_walls, int dim, int N_res, double* U) : variables(N, N_walls, dim, N_res)
{
    for(uint n = 0; n < N; n++)  
    {  
        for(uint i = 0; i < dim; i++)
        {
            W(n,i) = U[i];
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
        p[i] = thermo::pressure(par,W(i));
    }
}

void variables::temperature(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        T[i] = thermo::temperature(par,W(i));
    }
}

void variables::mach_number(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        M[i] = thermo::mach_number_stagnate(par,W(i));
    }
}
