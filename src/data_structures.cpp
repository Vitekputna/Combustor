#include "data_structures.h"
#include <stdlib.h>

variables::variables(int N, int dim) : N{N}, dim{dim}
{
    N_max = 6*N+dim*N;
    mem_ptr = (double*)(calloc(N_max,sizeof(double)));

    rho = mem_ptr;
    rhou = rho + N;
    rhov = rhou + N;
    e = rhov + N;
    p = e + N;
    T = p + N;
    wall_flux = T + N;
}

variables::variables(int N, int dim, double* U) : variables(N, dim)
{

    for(uint i = 0; i < N; i++)
    {
        rho[i] = U[0];
        rhou[i] = U[1];
        rhov[i] = U[2];
        e[i] = U[3];
    }
}

variables::~variables()
{
    free(mem_ptr);
}