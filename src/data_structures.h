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

<<<<<<< HEAD
    unsigned int iter,n_t = 50,n_r = 1000,n_b = 100;
=======
    unsigned int iter,n_t = 100,n_r = 1000,n_b = 100;
>>>>>>> 2aa42f5a4cf211ee82c62dc741fb2c058266c72a

    int res_idx = 3;
    double max_res = 1e-3; 
    int min_iter = 100;
    int max_iter = 1e6;
    int bisec_iter = 5;
};

struct variables
{
    int N, N_walls, dim;

    double *p, *T, *M, *res;

    array W,wall_flux;

    variables(int N, int N_walls, int dim);
    variables(int N, int N_walls, int dim, double* U);
    ~variables();

    void pressure(parameters const& par);
    void temperature(parameters const& par);
    void mach_number(parameters const& par);
};