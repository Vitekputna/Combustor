#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include "mesh.h"
#include "math.h"

typedef unsigned int uint;

face::face(vec1d const& a, vec1d const& b)
{
    double sx,sy;

    sx = b[0]-a[0];
    sy = b[1]-a[1];

    xf = a[0] + sx/2;
	yf = a[1] + sy/2; 

    S = sqrt(sx*sx + sy*sy);

    n[0] = sy/S;
    n[1] = -sx/S;
    s[0] = sx;
    s[1] = sy;
}

face::face(){}

cell::cell()
{
    x = 0, y = 0;
}

cell::cell(vec1ui nodes, vec2d const& all_nodes)
{
    std::vector<double> X,Y;
    X.reserve(4);
    Y.reserve(4);

    for(int i = 0; i < 4; i++)
    {
        X[i] = all_nodes[nodes[i]][0];
        Y[i] = all_nodes[nodes[i]][1];

        x += 0.25*X[i];
        y += 0.25*Y[i];
    }

    V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) )
      + 0.5*abs( (X[3] - X[2]) * (Y[0] - Y[3]) - (X[0] - X[3]) * (Y[3] - Y[2]) ); 
}

void cell::add_cell_wall(unsigned int wall_idx)
{
    cell_walls[free_wall_slot_idx] = wall_idx;
    free_wall_slot_idx++;
}

group::group()
{
    member_idx = std::vector<uint>{};
}

mesh::mesh(int i)
{
    name = "";
}

mesh::mesh(std::string path)
{
    name = path.substr(0,path.find('.'));
    std::cout << "Loading mesh: " << name << "\n";

    load_mesh(path,nodes,edges,quads);
    N_cells = quads.size();
    construct_ghost_cells();
    N_ghosts = quads.size() - N_cells;

    cells.resize(N_cells+N_ghosts);

    construct_cells();
    sort_mesh();
    set_owner_idx();
    N_walls = walls.size();
    N = quads.size();

    group_inlets();
    
    std::cout << "Mesh loaded, number of walls: " << N_walls << " , number of cells: " 
              << N_cells << " number of ghosts: " << N_ghosts << "\n\n";
}

void mesh::load_mesh(std::string path, vec2d& nodes, vec2ui& edges, vec2ui& quads)
{
    std::ifstream file(path);
    std::cout << "Warning: This reader takes in gmsh format, counting from 1 not zero!!!\n";

	if(file.fail())
	{
		std::cout << "//////////////////////////////////// \n";
		std::cout << "File not found \n";
		std::cout << "//////////////////////////////////// \n";
		throw std::exception();
	}

	std::vector<std::string> text_vec;
	std::string text;

	// Reading as text file
	while (std::getline(file, text))
	{
		text_vec.push_back(text);
	}

	int num_nodes = stoi(text_vec[4]);

	// Reading nodes
	std::string word = "";
	double number;
    int integer;
	char k;
	std::vector<double> row;
    std::vector<unsigned int> rowi;

	//std::cout << "NODES" << std::endl;
	for (unsigned int i = 5; i < num_nodes+5; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				number = (std::stod(word));
				row.push_back(number);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		nodes.push_back(row);
		row.clear();
		//std::cout << "\n";
	}

	// reading edge 
	//std::cout << "EDGES" << std::endl;
	int edge_start = num_nodes + 6;
	int num_edges = stoi(text_vec[edge_start]);

	//std::cout << num_edges << "\n";

	for (unsigned int i = edge_start+1; i < num_edges + edge_start+1; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				integer = std::stoi(word);
				rowi.push_back(integer-1);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		edges.push_back(rowi);
		rowi.clear();
		//std::cout << "\n";
	}

	// reading quads
	//std::cout << "QUADS" << std::endl;
	int quads_start = num_edges + edge_start + 2;
	int num_quads = stoi(text_vec[quads_start]);

	//std::cout << num_quads << "\n";

	for (unsigned int i = quads_start+1; i < quads_start + num_quads +1; i++)
	{
		for (unsigned int j = 0; j < text_vec[i].length(); j++)
		{
			if (text_vec[i][j] != 32)
			{
				while (text_vec[i][j] != 32 && j < text_vec[i].length())
				{
					k = text_vec[i][j];
					word.push_back(k);
					j++;
				}
				integer = std::stoi(word);
				rowi.push_back(integer-1);
				//std::cout << row.back() << " ";
				word = "";
			}
		}
		quads.push_back(rowi);
		rowi.clear();
		//std::cout << "\n";
	}
}

void mesh::print_mesh()
{
    for(auto const& node : nodes)
    {
        for(auto const& i : node)
        {
            std::cout << i << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";


    for(auto const& quad : quads)
    {
        for(auto const& i : quad)
        {
            std::cout << i << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    for(auto const& edge : edges)
    {
        for(auto const& i : edge)
        {
            std::cout << i << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int mesh::find_neigbour_cell(vec1i common_nodes_idx, int owner_idx)
{
    bool found_1 = false, found_2 = false;

    unsigned int quad;
    for(quad = 0; quad < quads.size(); quad++)
    {
        for(auto const& q : quads[quad])
        {
            if(q == common_nodes_idx[0])
            {
                found_1 = true;
            }
            else if (q == common_nodes_idx[1])
            {
                found_2 = true;
            }
        }

        if(found_1 && found_2 && quad != owner_idx){return quad;}

        found_1 = false;
        found_2 = false;
    }
    return -1;
}

bool mesh::wall_uniqueness(face const& new_wall)
{
    for(auto const& wall : walls)
    {
        if(wall.owner_cell_index == new_wall.owner_cell_index && wall.neigbour_cell_index == new_wall.neigbour_cell_index ||
           wall.owner_cell_index == new_wall.neigbour_cell_index && wall.neigbour_cell_index == new_wall.owner_cell_index)
        {
            return false;
        }
    }
    return true;
}

void mesh::construct_ghost_cells()
{
     vec2ui ghost_quads;
     quads_offset = quads.size();
    
     int i = 0;
     for(auto const& edge : edges)
     {  
        if(edge.back() != -1)
        {
            ghost_quads.push_back(std::vector<unsigned int>{edge[0],edge[1],edge[0],edge[1]});
            ghost_cell_idx.push_back(i+quads_offset);
            ghost_cell_val.push_back(edge.back());
            i++;
        }
         
     }

     quads.insert(quads.end(),ghost_quads.begin(),ghost_quads.end());
}

void mesh::sort_mesh()
{
    int neighbour;

    std::vector<std::vector<int>> node_idx = {{0,1},{1,2},{2,3},{3,0}};
    std::vector<int> node_vec;

    unsigned int wall_idx = 0;
    for(unsigned int c_idx = 0; c_idx < quads.size(); c_idx++)
    {
        for(unsigned int w = 0; w < 4; w++)
        {
            node_vec = {int(quads[c_idx][node_idx[w][0]]),int(quads[c_idx][node_idx[w][1]])};
            neighbour = find_neigbour_cell(node_vec,c_idx);

            face wall(nodes[node_vec[0]],nodes[node_vec[1]]);
            wall.owner_cell_index = c_idx;
            wall.neigbour_cell_index = neighbour;

            

            if(wall_uniqueness(wall) && neighbour != -1)
            {
                walls.push_back(wall);
                //cells[c_idx].owner_idx[cells[c_idx].free_wall_slot_idx] = 1;    
                cells[c_idx].add_cell_wall(wall_idx);
                
                cells[neighbour].add_cell_wall(wall_idx);
                wall_idx++;
            }
        }
    }
}

void mesh::set_owner_idx()
{
    int i = 0;
    for(auto& cell : cells)
    {
        int j = 0;
        for(auto const& w : cell.cell_walls)
        {
            if(walls[w].owner_cell_index == i)
            {
                cell.owner_idx[j] = 1;
            }
            j++;
        }
        i++;
    }
}

void mesh::construct_cells()
{
    for(uint k = 0; k < quads.size();k++)
    {
        //cells.push_back(cell(quads[k], nodes));
        cells[k] = cell(quads[k],nodes);

        if(cells[k].V > 0)
        {
            min_V = std::min(min_V,cells[k].V);
        }
        
    }
}

void mesh::group_inlets()
{
    std::vector<int> group_idx;
    std::vector<int> group_size;

    int max_group_size = 0;

    for(uint c = 0; c < ghost_cell_val.size(); c++)
    {
        auto val = ghost_cell_val[c];

        if(val != -1 && 
           !std::count(group_idx.begin(),group_idx.end(),val))
        {
            group_idx.push_back(val);
            //std::cout << "added: " << val << "\n";
            group_size.push_back(0);
            boundary_groups.push_back(group());
        }

        if(val != -1)
        {
            for(uint i = 0; i < group_idx.size(); i++)
            {
                if(val == group_idx[i])
                {
                    group_size[i]++;
                    boundary_groups[i].member_idx.push_back(ghost_cell_idx[c]);
                    boundary_groups[i].group_value = val;
                }
            }
        }
    }

    int g = 0;
    std::cout << "Boundary groups: \n";
    for(auto const& group : boundary_groups)
    {
        std::cout << "group: " << g << " ,value: " << group.group_value
                  << " ,size: " << group.member_idx.size() << "\n";
        g++;
    }

    // int g = 0;
    // for(auto const& group : boundary_groups) 
    // {
    //     std::cout << "Group: " << g << "\n"; 

    //     for(auto const& item : group.member_idx)
    //     {
    //         std::cout << item << " : " << group.group_value << "\n";
    //     }

    //     g++;
    // }
    // std::cout << "\n";
}

void mesh::export_mesh()
{
    std::ofstream f(name + "_walls.txt");
    f << N_walls << "\n";
    for(auto const& wall : walls)
    {
        f << wall.xf << " " << wall.yf << " " << wall.S << "\n";
        f << wall.n[0] << " " << wall.n[1] << "\n";
        f << wall.s[0] << " " << wall.s[1] << "\n";
        f << wall.owner_cell_index << " " << wall.neigbour_cell_index << "\n";
    }
    f.close();

    std::ofstream ff(name + "_cells.txt");
    ff << N_cells << " " << N_ghosts << "\n";
    for(auto const& cell : cells)
    {
        ff << cell.x << " " << cell.y << " " << cell.V << "\n";
        ff << cell.cell_walls[0] << " " << cell.cell_walls[1] 
           << " " << cell.cell_walls[2] << " " << cell.cell_walls[3] << "\n";
        ff << cell.owner_idx[0] << " " << cell.owner_idx[1] << " "
           << cell.owner_idx[2] << " " << cell.owner_idx[3] << "\n";
    }
    ff.close();
}

std::vector<double> mesh::extract(std::string& text)
{
    std::string word;
    std::vector<double> res;

    for (unsigned int j = 0; j < text.length(); j++)
    {
        if (text[j] != 32)
        {
            while (text[j] != 32 && j < text.length())
            {
                word.push_back(text[j]);
                j++;
            }
            res.push_back(std::stod(word));
            word = "";
        }
    }
    
    return res;
}

void mesh::import_mesh(std::string path)
{   
    name = path.substr(0,path.find('.'));
    std::cout << "Loading mesh: " << name << "\n";

    load_mesh(path,nodes,edges,quads);


    std::ifstream file("mesh/" + name +"_cells.txt");

    if(file.fail())
	{
		std::cout << "//////////////////////////////////// \n";
		std::cout << "File not found \n";
		std::cout << "//////////////////////////////////// \n";
		throw std::exception();
	}

    std::string text;

    std::getline(file,text);

    std::vector<std::string> text_vec;
    std::vector<double> double_vec;

    double_vec = extract(text);

    N_cells = uint(double_vec[0]);
    N_ghosts = uint(double_vec[1]);
    N = N_cells + N_ghosts;

    for(unsigned int i = 0; i < N_cells + N_ghosts; i++)
    {
        cell C;

        std::getline(file, text);
        double_vec = extract(text);
        C.x = double_vec[0];
        C.y = double_vec[1];
        C.V = double_vec[2];

        std::getline(file, text);
        double_vec = extract(text);
        C.cell_walls[0] = double_vec[0];
        C.cell_walls[1] = double_vec[1];
        C.cell_walls[2] = double_vec[2];
        C.cell_walls[3] = double_vec[3];

        std::getline(file, text);
        double_vec = extract(text);
        C.owner_idx[0] = double_vec[0];
        C.owner_idx[1] = double_vec[1];
        C.owner_idx[2] = double_vec[2];
        C.owner_idx[3] = double_vec[3];

        cells.push_back(C);
    }

    file.close();

    std::ifstream f("mesh/" + name +"_walls.txt");

    if(f.fail())
	{
		std::cout << "//////////////////////////////////// \n";
		std::cout << "File not found \n";
		std::cout << "//////////////////////////////////// \n";
		throw std::exception();
	}

    std::getline(f,text);

    double_vec = extract(text);

    N_walls = uint(double_vec[0]);

    for(unsigned int i = 0; i < N_walls; i++)
    {
        face C;

        std::getline(f, text);
        double_vec = extract(text);
        C.xf = double_vec[0];
        C.yf = double_vec[1];
        C.S  = double_vec[2];

        std::getline(f, text);
        double_vec = extract(text);
        C.n[0] = double_vec[0];
        C.n[1] = double_vec[1];

        std::getline(f, text);
        double_vec = extract(text);
        C.s[0] = double_vec[0];
        C.s[1] = double_vec[1];

        std::getline(f, text);
        double_vec = extract(text);
        C.owner_cell_index = double_vec[0];
        C.neigbour_cell_index = double_vec[1];

        walls.push_back(C);
    }

    file.close();


}