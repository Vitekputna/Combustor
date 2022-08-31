#include <fstream>
#include <string>
#include "data_structures.h"
#include "mesh.h"

typedef unsigned int uint;

void export_vtk(variables& var,mesh const& MESH, std::string name)
{
	std::ofstream f(name);

	std::cout << "Exporting to paraview format...\n";

	f << "# vtk DataFile Version 2.0" << std::endl;
	f << "Unstructured Grid Example" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET UNSTRUCTURED_GRID" << std::endl;
	f << std::endl;
	f << "POINTS" << " " << MESH.nodes.size() << " " << "float" << std::endl;

	for (unsigned int i = 0; i < MESH.nodes.size(); i++)
	{
		for (unsigned int j = 0; j < MESH.nodes[0].size(); j++)
		{
			f << MESH.nodes[i][j] << " ";
		}
		f << std::endl;
	}

	f << std::endl;
	f << "CELLS " << MESH.N_cells << " " << MESH.N_quads * 5 + MESH.N_trigs * 4 << std::endl;
	
    for (unsigned int i = 0; i < MESH.N_quads; i++)	
    {
        f << "4 ";
        for(unsigned int k = 0; k < MESH.quads[0].size();k++)
        {
            f << MESH.quads[i][k] << " ";
        }
        f << "\n";
    }
	for (unsigned int i = 0; i < MESH.N_trigs; i++)	
    {
        f << "3 ";
        for(unsigned int k = 0; k < MESH.trigs[0].size();k++)
        {
            f << MESH.trigs[i][k] << " ";
        }
        f << "\n";
    }


	f << std::endl;
	f << "CELL_TYPES " << MESH.N_cells << std::endl;
	for (unsigned int i = 0; i < MESH.N_quads; i++)
	{
		f << "9" << std::endl;
	}
	for (unsigned int i = 0; i < MESH.N_trigs; i++)
	{
		f << "5" << std::endl;
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

	std::cout << "Done! \n\n";
}

void export_res(variables& var, std::string name)
{
	std::ofstream f(name);

	for(uint i = 0; i < (uint)(1e6/3000+100); i++)
	{
		f << var.res[i] << "\n";
	}
}