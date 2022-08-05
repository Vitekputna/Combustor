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
    double* operator()(int i);
};

struct parameters
{
    double gamma = 1.4;
    double r = 287;
};

struct config
{
    double dt = 1e-50;

    double CFL = 1;
    unsigned int iter,n_t;
    
};

struct variables
{
    int N, N_walls, dim;

    double *p, *T;

    array W,wall_flux;

    variables(int N, int N_walls, int dim);
    variables(int N, int N_walls, int dim, double* U);
    ~variables();

    void pressure(parameters const& par);
    void temperature(parameters const& par);
};