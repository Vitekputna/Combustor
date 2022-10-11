#pragma once
#include "data_structures.h"
#include "mesh.h"

void reconstruct(int dim, int c, int w, double* Wr, double* W, double* grad, mesh const& msh)
{
    double dx = msh.walls[c].xf - msh.cells[w].x;
    double dy = msh.walls[c].yf - msh.cells[w].y;

    for(int i = 0; i < dim; i++)
    {
        Wr[i] = W[i];// + grad[2*i]*dx + grad[2*i+1]*dy;
    }
}

void compute_wall_flux(variables& var, mesh const& msh, parameters const& par, void(*flux)(int dim,
                                                                                           double* w,
                                                                                           double* n,
                                                                                           double* o,
                                                                                           parameters const&,
                                                                                           face const&))
{
    int n,o;

    double Wo[var.dim];
    double Wn[var.dim];

    #pragma omp parallel for private(o,n,Wn,Wo) shared(msh,var,par)
        for(int w = 0; w < msh.N_walls;w++)
        {
            n = msh.walls[w].neigbour_cell_index;
            o = msh.walls[w].owner_cell_index;

            reconstruct(var.dim,o,w,Wo,var.W(o),var.grad(o),msh);
            reconstruct(var.dim,n,w,Wn,var.W(n),var.grad(n),msh);

            flux(var.dim,var.wall_flux(w),Wn,Wo,par,msh.walls[w]);
        }
}

void compute_cell_res(variables& var, mesh const& msh, config const& cfg,void(*source_func)(double,double*))
{
    int f;
    #pragma omp parallel for private(f) shared(var,cfg,msh)
    for(int c = 0; c < msh.N_cells; c++)
    {
        for(uint k = 0; k < var.dim; k++)
        {
            f = 0;
            for(auto const& wall : msh.cells[c].cell_walls)
            {
                var.W(c,k) += -cfg.dt/msh.cells[c].V*var.wall_flux(wall,k)*msh.cells[c].owner_idx[f];
                f++;
            }
            
        }
        source_func(cfg.dt,var.W(c));
    }
}

void compute_cell_gradient(variables& var, mesh const& msh)
{
    int o,n;
    double V;
    #pragma omp parallel for private(o,n,V) shared(var,msh)
    for(int c = 0; c < msh.N_cells;c++)
    {   
        V = msh.cells[c].V;

        for(int k = 0; k < var.dim; k++) // Neslo by to prehodit s loopem pres steny?
        {   
            int idx = 0;
            for(auto const& w : msh.cells[c].cell_walls)
            {
                n = msh.walls[w].neigbour_cell_index;
                o = msh.walls[w].owner_cell_index;

                o = (o == c) ? n : o;

                var.grad(c,2*k) += 1/V*var.W(o,k)*msh.walls[w].n[0]*msh.walls[w].S*msh.cells[c].owner_idx[idx];
                var.grad(c,2*k+1) += 1/V*var.W(o,k)*msh.walls[w].n[1]*msh.walls[w].S*msh.cells[c].owner_idx[idx];

                idx++;
            }
        }
    }
}

void grad_limiting(variables& var, mesh const& msh)
{
    int o,n;
    double alfa = 1;
    double d_max, d_min;
    for(int c = 0; c < msh.N_cells; c++)
    {
        for(int k = 0; k < var.dim; k++)
        {
            d_max = -INFINITY, d_min = INFINITY;
            for(auto const& w : msh.cells[c].cell_walls)
            {
                n = msh.walls[w].neigbour_cell_index;
                o = msh.walls[w].owner_cell_index;
                n = (n == c) ? o : n;

                double dx = msh.walls[c].xf - msh.cells[w].x;
                double dy = msh.walls[c].yf - msh.cells[w].y;

                d_max = std::max(d_max,var.grad(c,2*k)*dx + var.grad(c,2*k+1)*dy);
                d_min = std::min(d_min,var.grad(c,2*k)*dx + var.grad(c,2*k+1)*dy);
            }


            alfa = 1;
            for(auto const& w : msh.cells[c].cell_walls)
            {
                n = msh.walls[w].neigbour_cell_index;
                o = msh.walls[w].owner_cell_index;
                n = (n == c) ? o : n;

                alfa = std::min(alfa, std::min(abs((var.W(n,k)-var.W(c,k))/(d_max)),
                                               abs((var.W(n,k)-var.W(c,k))/(d_min))));
            }

            var.grad(c,2*k) *= 0.5*alfa;
            var.grad(c,2*k+1) *= 0.5*alfa;
            var.alfa(c,k) = 0.5*alfa;
        }
    }
}