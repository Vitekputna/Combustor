#pragma once
#include "data_structures.h"
#include "mesh.h"
#include "math.h"
#include "omp.h"

extern int N_threads;

void reconstruct(int dim, int c, int w, double* Wr, double* W, double* grad, mesh const& msh)
{
    double dx = msh.walls[c].xf - msh.cells[w].x;
    double dy = msh.walls[c].yf - msh.cells[w].y;

    for(int i = 0; i < dim; i++)
    {
        Wr[i] = W[i];// + grad[2*i]*dx + grad[2*i+1]*dy;
    }
}

void compute_wall_flux(variables& var, mesh const& msh, std::vector<parameters> const& par, void(*flux)(int,int,
                                                                                           double*,
                                                                                           double*,
                                                                                           double*,
                                                                                           std::vector<parameters> const&,
                                                                                           face const&))
{
    int n,o;
    double Wo[var.dim];
    double Wn[var.dim];

    #pragma omp parallel num_threads(N_threads) private(o,n,Wn,Wo) shared(msh,var,par)
    {
        #pragma omp for
        for(int w = 0; w < msh.N_walls;w++)
        {
            n = msh.walls[w].neigbour_cell_index;
            o = msh.walls[w].owner_cell_index;

            reconstruct(var.dim,o,w,Wo,var.W(o),var.grad(o),msh);
            reconstruct(var.dim,n,w,Wn,var.W(n),var.grad(n),msh);

            flux(var.vel_comp,var.n_comp,var.wall_flux(w),Wn,Wo,par,msh.walls[w]);
        }
    }
}

void compute_diffusive_flux(variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    for(int i = 0; i < msh.N_walls; i++)
    {
        var.diff_flux(i,4) = par[0].lambda*msh.walls[i].S*var.T_grad[i];
    }
}

void apply_cell_res(variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par,
                      void(*source_func)(array&,variables&,mesh const&,config const&,std::vector<parameters> const&))
{
    int f;
    #pragma omp parallel private(f) shared(var,cfg,msh)
    {
        #pragma omp for
        for(int c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                f = 0;
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    var.W(c,k) += -cfg.dt/msh.cells[c].V*(var.wall_flux(wall,k)*msh.cells[c].owner_idx[f]);
                    f++;
                }
            }
        }

        source_func(var.W,var,msh,cfg,par);
    }
}

void compute_cell_res(array& res ,variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par,
                      void(*source_func)(array&,variables&,mesh const&,config const&,std::vector<parameters> const&))
{
    int f;
    #pragma omp parallel private(f) shared(var,cfg,msh)
    {
        #pragma omp for
        for(int c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                f = 0;
                for(auto const& wall : msh.cells[c].cell_walls)
                {
                    res(c,k) += -1/msh.cells[c].V*(var.wall_flux(wall,k)*msh.cells[c].owner_idx[f]);
                    f++;
                }
            }
        }
        source_func(res,var,msh,cfg,par);
    }
}

void apply_cell_res(array& old_res, array& new_res, double old_coef, double new_coef, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    int f;
    #pragma omp parallel private(f) shared(var,cfg,msh)
    {
        #pragma omp for
        for(int c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                var.W(c,k) += cfg.dt*old_coef*old_res(c,k) + cfg.dt*new_coef*new_res(c,k);
                old_res(c,k) = 0;
                new_res(c,k) = 0;
            }
        }
    }
}

void apply_cell_res(array& res, variables& var, mesh const& msh, config const& cfg, std::vector<parameters> const& par)
{
    int f;
    #pragma omp parallel private(f) shared(var,cfg,msh)
    {
        #pragma omp for
        for(int c = 0; c < msh.N_cells; c++)
        {
            for(uint k = 0; k < var.dim; k++)
            {
                var.W(c,k) += cfg.dt*res(c,k);
                res(c,k) = 0;
            }
        }
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

        for(auto const& w : msh.cells[c].cell_walls)
        {
            int idx = 0;

            for(int k = 0; k < var.dim; k++)
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

void compute_wall_gradiend(int j, variables& var, mesh const& msh)
{
    double Nx,Ny;
    double nx,ny;
    double fx,fy;
    double xf,yf;
    double S,V_cell;

    std::vector<std::vector<double>> V(5,std::vector<double>(2,0));
    std::vector<double> U(5,0.0);

    int w = 0;
    #pragma omp parallel for private(V,U,Nx,Ny,nx,ny,xf,yf,S,V_cell) shared(msh,var)
    for(auto const& wall : msh.walls)
    {
        xf = wall.xf;
        yf = wall.yf;

        V[0][0] = msh.cells[wall.owner_cell_index].x;
        V[0][1] = msh.cells[wall.owner_cell_index].y;

        V[1][0] = msh.nodes[wall.vertices[0].node_idx][0];
        V[1][1] = msh.nodes[wall.vertices[0].node_idx][1];

        V[2][0] = msh.cells[wall.neigbour_cell_index].x;
        V[2][1] = msh.cells[wall.neigbour_cell_index].y;

        V[3][0] = msh.nodes[wall.vertices[1].node_idx][0];
        V[3][1] = msh.nodes[wall.vertices[1].node_idx][1];

        V[4][0] = msh.cells[wall.owner_cell_index].x;
        V[4][1] = msh.cells[wall.owner_cell_index].y;

        U[0] = var.W(wall.owner_cell_index,j);
        U[1] = var.compute_vertex_average(wall.vertices[0],msh)[j];
        U[2] = var.W(wall.neigbour_cell_index,j);
        U[3] = var.compute_vertex_average(wall.vertices[1],msh)[j];
        U[4] = var.W(wall.owner_cell_index,j);

        V_cell = 0.5*(V[0][1]*V[1][0] - V[0][0]*V[1][1]) + 0.5*(V[2][1]*V[3][0] - V[2][0]*V[3][1]);

        for(int i = 0; i < 3; i++)
        {
            Nx = V[i+1][1] - V[i][1];
            Ny = V[i][0] - V[i+1][0];

            S = sqrt(Nx*Nx+Ny*Ny);

            fx = 0.5*(V[i+1][0] + V[i][0]) - xf;
            fy = 0.5*(V[i+1][1] + V[i][1]) - yf;

            nx = Nx/(S)*(Nx*fx+Ny*fy)/(abs(Nx*fx+Ny*fy));
            ny = Ny/(S)*(Nx*fx+Ny*fy)/(abs(Nx*fx+Ny*fy));

            var.wall_grad(w,j) += 1/V_cell*(nx*wall.n[0] + ny*wall.n[1])*S*(U[i+1] - U[i])/2;
        }

        w++;
    }
}

void compute_wall_T_gradiend(variables& var, mesh const& msh)
{
    double Nx,Ny;
    double nx,ny;
    double fx,fy;
    double xf,yf;
    double S,V_cell;

    std::vector<std::vector<double>> V(5,std::vector<double>(2,0));
    std::vector<double> U(5,0.0);

    int w = 0;
    for(auto const& wall : msh.walls)
    {
        xf = wall.xf;
        yf = wall.yf;

        V[0][0] = msh.cells[wall.owner_cell_index].x;
        V[0][1] = msh.cells[wall.owner_cell_index].y;

        V[1][0] = msh.nodes[wall.vertices[0].node_idx][0];
        V[1][1] = msh.nodes[wall.vertices[0].node_idx][1];

        V[2][0] = msh.cells[wall.neigbour_cell_index].x;
        V[2][1] = msh.cells[wall.neigbour_cell_index].y;

        V[3][0] = msh.nodes[wall.vertices[1].node_idx][0];
        V[3][1] = msh.nodes[wall.vertices[1].node_idx][1];

        V[4][0] = msh.cells[wall.owner_cell_index].x;
        V[4][1] = msh.cells[wall.owner_cell_index].y;

        U[0] = var.T[wall.owner_cell_index];
        U[1] = var.compute_T_vertex_average(wall.vertices[0],msh);
        U[2] = var.T[wall.neigbour_cell_index];
        U[3] = var.compute_T_vertex_average(wall.vertices[1],msh);
        U[4] = var.T[wall.owner_cell_index];

        V_cell = 0.5*(V[0][1]*V[1][0] - V[0][0]*V[1][1]) + 0.5*(V[2][1]*V[3][0] - V[2][0]*V[3][1]);

        for(int i = 0; i < 3; i++)
        {
            Nx = V[i+1][1] - V[i][1];
            Ny = V[i][0] - V[i+1][0];

            S = sqrt(Nx*Nx+Ny*Ny);

            fx = 0.5*(V[i+1][0] + V[i][0]) - xf;
            fy = 0.5*(V[i+1][1] + V[i][1]) - yf;

            nx = Nx/(S)*(Nx*fx+Ny*fy)/(abs(Nx*fx+Ny*fy));
            ny = Ny/(S)*(Nx*fx+Ny*fy)/(abs(Nx*fx+Ny*fy));

            var.T_grad[w] += 1/V_cell*(nx*wall.n[0] + ny*wall.n[1])*S*(U[i+1] - U[i])/2;
        }

        w++;
    }
}