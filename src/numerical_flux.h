#pragma once
#include "data_structures.h"
#include "mesh.h"

void HLL_flux(int vel_comp, int n_comp, double* w, double* n, double* o, parameters const& par, face const& f);
void HLL_flux_axi(int vel_comp, int n_comp, double* w, double* n, double* o, parameters const& par, face const& f);