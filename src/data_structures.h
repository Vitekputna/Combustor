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
    int N, N_walls, dim;

    double *p, *T;

    array W,wall_flux;

    variables(int N, int N_walls, int dim);
    variables(int N, int N_walls, int dim, double* U);
    ~variables();
};

struct parameters
{
    double gamma = 1.4;
    double r = 287;
};

struct config
{
    double dt, CFL;
    unsigned int iter,n_t;
    
};