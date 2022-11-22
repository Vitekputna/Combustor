#pragma once
#include "data_structures.h"
#include "mesh.h"

void HLL_flux(int vel_comp, int n_comp, double* w, double* n, double* o, std::vector<parameters> const& par, face const& f);