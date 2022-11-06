#pragma once
#include "mesh.h"
#include <vector>

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
    //dimenze ulohy
    int dim = 4;
    int vel_comp = 2;

    //double dt = 1e-50;
    double dt = 1e-7;
    double CFL = 1;

    unsigned int iter,n_t = 50,n_r = 1000,n_b = 100, n_exp = 100000;
    int res_idx = 3;
    double max_res = 1e-3;
    int min_iter = 100;
    int max_iter = 3e6;
    int bisec_iter = 5;

    double export_interval = 0.1;
    double max_time = 1;

    //axisymetric stuff
    double b = 5;
    int r_variable_idx = 1;

};

struct variables
{
    int N, N_walls, N_res;
    int dim, vel_comp;

    double *p, *T, *M, *res;

    array W,wall_flux,grad,alfa,Source;

    variables(int N, int N_walls, int dim, int vel_comp, int N_res);
    variables(int N, int N_walls, int dim, int vel_comp, int N_res, std::vector<double>& U);
    variables(mesh const msh, config const cfg, std::vector<std::vector<double>> U);
    ~variables();

    void pressure(parameters const& par);
    void temperature(parameters const& par);
    void mach_number(parameters const& par);

    std::vector<double> compute_vertex_average(vertex const& node, mesh const& msh);
};