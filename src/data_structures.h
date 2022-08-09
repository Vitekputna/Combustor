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

    unsigned int iter,n_t,n_r;

    int res_idx = 3;
    double max_res = 1e-3; 
    int min_iter = 1e3;
    int max_iter = 1e9;
};

struct variables
{
    int N, N_walls, dim;

    double *p, *T, *M;

    array W,wall_flux;

    variables(int N, int N_walls, int dim);
    variables(int N, int N_walls, int dim, double* U);
    ~variables();

    void pressure(parameters const& par);
    void temperature(parameters const& par);
    void mach_number(parameters const& par);
};