#include "data_structures.h"
#include <stdlib.h>

variables::variables(int N, int dim) : N{N}, dim{dim}
{
    N_max = 6*N+dim*N;
    mem_ptr = (double*)(calloc(N_max,sizeof(double)));

    rho = mem_ptr;
    u = rho + N;
    v = u + N;
    e = v + N;
    p = e + N;
    T = p + N;
    wall_flux = T + N;
}

variables::~variables()
{
    free(mem_ptr);
}