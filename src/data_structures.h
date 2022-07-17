#pragma once

struct array
{
    int k;
    double* arr;
    array();
    void allocate(int N, int m_k);
    array(int N,int k);
    ~array();
    double& operator()(int i,int j);
};

struct variables
{
    int N, dim, N_max;

    double *p, *T;

    //double* mem_ptr;
    //double* wall_flux_ptr;

    array W,wall_flux;

    variables(int N, int N_walls, int dim);
    variables(int N, int N_walls, int dim, double* U);
    ~variables();
};

struct parameters
{
    double gamma, r;
};

struct config
{
    double dt, CFL;
};