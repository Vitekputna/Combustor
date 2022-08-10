#include <fstream>
#include <string>
#include "data_structures.h"
#include "mesh.h"

void export_vtk(variables& var,mesh const& MESH, std::string name)
{
	std::ofstream f(name);

	f << "# vtk DataFile Version 2.0" << std::endl;
	f << "Unstructured Grid Example" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	f << std::endl;
	f << "POINTS" << " " << MESH.nodes.size() << " " << "float" << std::endl;

	for (unsigned int i = 0; i < MESH.nodes.size(); i++)
	{
		for (unsigned int j = 0; j < MESH.nodes[0].size()-1; j++)
		{
			f << MESH.nodes[i][j] << " ";
		}
		f << std::endl;
	}

	f << std::endl;
	f << "CELLS " << MESH.quads.size()-MESH.N_ghosts << " " << (MESH.quads.size()-MESH.N_ghosts) * 5 << std::endl;
	
    for (unsigned int i = 0; i < MESH.quads.size()-MESH.N_ghosts; i++)	
    {
        f << "4 ";
        for(unsigned int k = 0; k < MESH.quads[0].size()-1;k++)
        {
            f << MESH.quads[i][k] << " ";
        }
        f << "\n";
    }

	f << std::endl;
	f << "CELL_TYPES " << MESH.quads.size()-MESH.N_ghosts << std::endl;
	for (unsigned int i = 0; i < MESH.quads.size()-MESH.N_ghosts; i++)
	{
		f << "9" << std::endl;
	}
	f << std::endl;

	f << "CELL_DATA " << MESH.N_cells << std::endl;
	f << "SCALARS " << "rho" << " float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)
	{
		f << var.W(j,0) << std::endl;
	}
	f << std::endl;

	f << "SCALARS " << "e" << " float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)	
	{
        f << var.W(j,3) << std::endl;
	}
	f << std::endl;

	f << "SCALARS " << "p" << " float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)
	{
		f << var.p[j] << std::endl;
	}
	f << std::endl;

	f << "SCALARS " << "M" << " float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)
    {
		f << var.M[j] << "\n";
	}
	f << std::endl;
 
	f << "SCALARS " << "T" << " float 1" << std::endl;
	f << "LOOKUP_TABLE default" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)
	{
		f << var.T[j] << std::endl;
	}
	f << std::endl;

	f << "VECTORS " << "velocity" << " float" << std::endl;
	for (unsigned int j = 0; j < MESH.N_cells; j++)
	{
		f << var.W(j,1) << " " << var.W(j,2) << " 0" << std::endl;
	}
	f << std::endl;
}