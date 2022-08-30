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

cell::cell(int N_walls, vec1ui nodes, vec2d const& all_nodes)
{
    N_faces = nodes.size();

    std::vector<double> X,Y;
    X.reserve(N_walls);
    Y.reserve(N_walls);

    for(int i = 0; i < N_walls; i++)
    {
        X[i] = all_nodes[nodes[i]][0];
        Y[i] = all_nodes[nodes[i]][1];

        x += 1/N_walls*X[i];
        y += 1/N_walls*Y[i];
    }

    if(N_walls == 4)
    {
        V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) )
            + 0.5*abs( (X[3] - X[2]) * (Y[0] - Y[3]) - (X[0] - X[3]) * (Y[3] - Y[2]) ); 
    }
    else
    {
        V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) );
    }
    
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
    N_quads = quads.size();
    N_trigs = trigs.size();

    construct_ghost_cells();
    N_ghosts = ghosts.size();

    N_cells = N_quads+N_trigs;

    cells.resize(N_cells+N_ghosts);

    construct_cells();
    sort_mesh();
    // set_owner_idx();
    // N_walls = walls.size();
    // N = quads.size();

    // group_inlets();
    
    // std::cout << "Mesh loaded, number of walls: " << N_walls << " , number of cells: " 
    //           << N_cells << " number of ghosts: " << N_ghosts << "\n\n";
}

template<typename T>
std::vector<std::vector<T>> mesh::read_segment(std::vector<std::string>& text,int from,int offset, int length)
{
    char k;
    double number;
    std::string word = "";
    std::vector<T> row;
    std::vector<std::vector<T>> res;

    int N = std::stoi(text[from]);
    from++;

    for (unsigned int i = from; i < from + N; i++)
	{
		//for (unsigned int j = 0; j < text[i].length(); j++)
        int n = 0;
        unsigned int j = 0;
        while(j < text[i].length() && n < length)
		{
			if (text[i][j] != 32)
			{
				while (text[i][j] != 32 && j < text[i].length())
				{
					k = text[i][j];
					word.push_back(k);
					j++;
				}

                n++;
				number = (T)(std::stod(word));
				row.push_back(number+offset);
				std::cout << row.back() << " ";
				word = "";
			}
        j++;
		}
		res.push_back(row);
		row.clear();
		std::cout << "\n";
	}

    return res;
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

	// Reading as text file na picu predelat na cteni a formatovani primo ze souboru
	while (std::getline(file, text))
	{
		text_vec.push_back(text);
	}

    std::string s1 = "Vertices";
    std::string s2 = "Edges";
    std::string s3 = "Quadrilaterals";
    std::string s4 = "Triangles";

    for(int i = 0; i < text_vec.size(); i++)
    {
        if(text_vec[i].find(s1) != std::string::npos)
        {
            std::cout << "reading Vertices...\n";
            nodes = read_segment<double>(text_vec,i+1,0,3);
            std::cout << nodes.size() << "\n";
        }
        else if(text_vec[i].find(s2) != std::string::npos)
        {
            std::cout << "reading Edges...\n";
            edges = read_segment<uint>(text_vec,i+1,-1,3);
            std::cout << edges.size() << "\n";
        }
        else if(text_vec[i].find(s3) != std::string::npos)
        {
            std::cout << "reading Quadrilaterals...\n";
            quads = read_segment<uint>(text_vec,i+1,-1,4);
            std::cout << quads.size() << "\n";
        }
        else if(text_vec[i].find(s4) != std::string::npos)
        {
            std::cout << "reading Triangles...\n";
            trigs = read_segment<uint>(text_vec,i+1,-1,3);
            std::cout << trigs.size() << "\n";
        }
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

    // Přidat mapu pro odlišení buněk které už mají všechny sousedy nalezené... možná optimalizace

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
    std::cout << "constructing ghosts...\n";

    N_ghosts = edges.size();
    ghosts.resize(N_ghosts,std::vector<uint>(4));
    ghost_cell_val.resize(N_ghosts);

    int i = 0;
    for(auto const& edge : edges)
    {  
       if(edge.back() != std::numeric_limits<uint>::max())
       {
            ghost_cell_val[i] = edge.back();
            ghosts[i] = (std::vector<unsigned int>{edge[0],edge[1],edge[0],edge[1]});
            i++;     
       }
    }


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
    for(uint k = 0; k < N_quads;k++)
    {
        cells[k] = cell(4,quads[k],nodes);

        if(cells[k].V > 0)
        {
            min_V = std::min(min_V,cells[k].V);
        }
    }

    for(uint k = 0; k < N_trigs;k++)
    {
        cells[k+N_quads] = cell(3,trigs[k],nodes);

        if(cells[k].V > 0)
        {
            min_V = std::min(min_V,cells[k].V);
        }
    }

    for(uint k = 0; k < N_ghosts;k++)
    {
        cells[k+N_quads+N_trigs] = cell(4,ghosts[k],nodes);
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
                    //boundary_groups[i].member_idx.push_back(ghost_cell_idx[c]);
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