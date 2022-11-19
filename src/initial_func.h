#pragma once
#include "data_structures.h"
#include "initial_cond.h"
#include "boundary.h"
#include "cmath"

void no_move_flow(int dim, initial_conditions& IC);
std::vector<double> move_flow(int vel_comp, int n_comp, initial_conditions& IC, std::vector<parameters> const& par);
void supersonic_inlet(int vel_comp, int n_comp, std::vector<parameters> const& par, boundary_group& bdr);
void supersonic_outlet(std::vector<parameters> const& par, boundary_group& bdr);
void subsonic_inlet(std::vector<parameters> const& par, boundary_group& bdr);
void subsonic_outlet(std::vector<parameters> const& par, boundary_group& bdr);