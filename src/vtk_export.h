#pragma once
#include "data_structures.h"
#include "mesh.h"
#include <string>

void export_vtk(variables& var,mesh const& MESH, std::string name);
void export_res(variables& var, std::string name);