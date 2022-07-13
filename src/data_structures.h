#pragma once

struct variables
{
    int N, dim, N_max;
    //double *rho, *u, *v, *e, *p, *T;
    double *rho, *rhou, *rhov, *e, *p, *T;
    double *wall_flux;

    double* mem_ptr;

    variables(int N, int dim);
    variables(int N, int dim, double* U);
    ~variables();
};

struct parameters
{
    double gamma, r;
};

struct config
{
    double dt, CFL;
    int dim;
};


