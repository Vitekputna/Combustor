#include "data_structures.h"
#include "thermodynamics.h"
#include <stdlib.h>

array::array(){}

array::array(int N, int k) : k{k}
{
    allocate(N,k);
}

void array::allocate(int N, int m_k)
{
    k = m_k;
    arr = (double*)(malloc(N*sizeof(double)));
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

variables::variables(int N, int N_walls, int dim) : N{N}, dim{dim}, N_walls{N_walls}
{
    W.allocate(N*dim,dim);
    wall_flux.allocate(N_walls*dim,dim);
}

variables::variables(int N, int N_walls, int dim, double* U) : variables(N, N_walls, dim)
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
}

void variables::pressure(parameters const& par)
{
    for(uint i = 0; i < N; i++)
    {
        //p[i] = pressure(par,W(i));
    }
}