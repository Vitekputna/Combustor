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
    vertices.reserve(2);

    double sx,sy;

    sx = b[0]-a[0];
    sy = b[1]-a[1];

    s[0] = sx;
    s[1] = sy;

    xf = a[0] + sx/2;
	yf = a[1] + sy/2; 

    S = sqrt(sx*sx + sy*sy);

    n[0] = sy/S;
    n[1] = -sx/S;
}

face::face(vec1d const& a, vec1d const& b, unsigned int ca, unsigned int cb, unsigned int ia, unsigned int ib)
{
    vertices.resize(2);

    vertices[0].node_idx = ia;
    vertices[0].common_cell_idx[0] = ca;
    vertices[0].common_cell_idx[1] = cb;

    vertices[1].node_idx = ib;
    vertices[1].common_cell_idx[0] = ca;
    vertices[1].common_cell_idx[1] = cb;

    double sx,sy;

    sx = b[0]-a[0];
    sy = b[1]-a[1];

    s[0] = sx;
    s[1] = sy;

    xf = a[0] + sx/2;
	yf = a[1] + sy/2; 

    S = sqrt(sx*sx + sy*sy);

    n[0] = sy/S;
    n[1] = -sx/S;
}

face::face(){}

cell::cell(){x = 0, y = 0;}

cell::cell(char N_walls, vec1ui nodes, vec2d const& all_nodes) : N_faces{N_walls}
{
    this->cell_walls.resize(N_walls);
    
    X.reserve(std::max((int)N_faces,3));
    Y.reserve(std::max((int)N_faces,3));

    for(int i = 0; i < N_walls; i++)
    {
        X[i] = all_nodes[nodes[i]][0];
        Y[i] = all_nodes[nodes[i]][1];

        x += X[i]/N_walls;
        y += Y[i]/N_walls;
    }

    if(N_faces == 4)
    {
        V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) )
            + 0.5*abs( (X[3] - X[2]) * (Y[0] - Y[3]) - (X[0] - X[3]) * (Y[3] - Y[2]) ); 
    }
    else if(N_faces == 3)
    {
        V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) );
    }
    else
    {
        V = 0;
    }
}

void cell::ghost(vec1ui nodes,vec2d const& all_nodes)
{
    for(int i = 0; i < 3; i++)
    {
        X[i] = all_nodes[nodes[i]][0];
        Y[i] = all_nodes[nodes[i]][1];

        x += X[i]/3;
        y += Y[i]/3;
    }
    V = 0.5*abs( (X[1] - X[0]) * (Y[2] - Y[1]) - (Y[1] - Y[0]) * (X[2] - X[1]) );
}

void cell::add_cell_wall(unsigned int wall_idx)
{
    cell_walls[free_wall_slot_idx] = wall_idx;
    free_wall_slot_idx++;
}

group::group(){member_idx = std::vector<uint>{};}

mesh::mesh(int i){name = "";}

mesh::mesh(std::string path)
{
    name = path.substr(0,path.find('.'));
    std::cout << "Loading mesh: " << name << "\n";

    load_mesh(path,nodes,edges,quads); //načtou se data ze souboru

    N_quads = quads.size(); // počet 4uhelníkových buněk
    N_trigs = trigs.size(); // počet 3uhelníkových buněk
    N_ghosts = edges.size(); // počet virtuálních buněk

    construct_ghost_cells(); // vytvoření dat k fiktivním buňkám

    N_cells = N_quads+N_trigs; // celokový počet reálných buněk
    N = N_quads+N_trigs+N_ghosts; // celkový počet buněk

    cells.resize(N_cells+N_ghosts);

    construct_cells();
    sort_mesh();
    set_owner_idx();
    N_walls = walls.size();
    
    group_inlets();

    expand_ghost_cells();

    vertex_init();
    
    std::cout << "Mesh loaded, number of walls: " << N_walls << " , number of cells: " 
              << N_cells << " number of ghosts: " << N_ghosts << "\n\n";
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
				//std::cout << row.back() << " ";
				word = "";
			}
        j++;
		}
		res.push_back(row);
		row.clear();
		//std::cout << "\n";
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

    std::vector<uint> quad_ps;
    std::vector<uint> trig_ps;

    std::vector<uint> ps_idx = {};

    for(int i = 0; i < text_vec.size(); i++)
    {
        if(text_vec[i].find(s1) != std::string::npos)
        {
            //std::cout << "reading Vertices...\n";
            nodes = read_segment<double>(text_vec,i+1,0,3);
            //std::cout << nodes.size() << "\n";
        }
        else if(text_vec[i].find(s2) != std::string::npos)
        {
            //std::cout << "reading Edges...\n";
            edges = read_segment<uint>(text_vec,i+1,-1,3);
            //std::cout << edges.size() << "\n";
        }
        else if(text_vec[i].find(s3) != std::string::npos)
        {
            //std::cout << "reading Quadrilaterals...\n";
            quads = read_segment<uint>(text_vec,i+1,-1,4);
            quad_ps.resize(quads.size());

            auto temp = read_segment<uint>(text_vec,i+1,0,5);
            for(int i = 0; i < temp.size(); i++)
            {
                quad_ps[i] = temp[i].back();
            }

            //std::cout << quads.size() << "\n";
        }
        else if(text_vec[i].find(s4) != std::string::npos)
        {
            //std::cout << "reading Triangles...\n";
            trigs = read_segment<uint>(text_vec,i+1,-1,3);
            trig_ps.resize(trigs.size());

            auto temp = read_segment<uint>(text_vec,i+1,0,4);
            for(int i = 0; i < temp.size(); i++)
            {
                trig_ps[i] = temp[i].back();
            }

            //std::cout << trigs.size() << "\n";
        }
    }

    ps_idx.resize(quad_ps.size() + trig_ps.size());

    for(int i = 0; i < quad_ps.size(); i++)
    {
        ps_idx[i] = quad_ps[i];
    }
    for(int i = quad_ps.size(); i < quad_ps.size()+trig_ps.size();i++)
    {
        ps_idx[i] = trig_ps[i];
    }

    group_surfaces(ps_idx);
}

void mesh::load_msh(std::string path, vec2d& nodes, vec2ui& edges, vec2ui& quads)
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

int mesh::find_neigbour_cell(vec2ui const& polygons, vec1i common_nodes_idx, int owner_idx)
{
    bool found_1 = false, found_2 = false;

    // Přidat mapu pro odlišení buněk které už mají všechny sousedy nalezené... možná optimalizace

    unsigned int pgon;
    for(pgon = 0; pgon < polygons.size(); pgon++)
    {
        for(auto const& p : polygons[pgon])
        {
            if(p == common_nodes_idx[0])
            {
                found_1 = true;
            }
            else if (p == common_nodes_idx[1])
            {
                found_2 = true;
            }
        }

        if(found_1 && found_2 && pgon != owner_idx){return pgon;}

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
    for(auto const& edge : edges)
    {
        if(edge.back() == std::numeric_limits<uint>::max())
        {
            N_ghosts--;
        }
    }

    ghosts.resize(N_ghosts,std::vector<uint>(2));
    ghost_cell_val.resize(N_ghosts);

    int i = 0;
    for(auto const& edge : edges)
    {  
       if(edge.back() != std::numeric_limits<uint>::max())
       {
            ghost_cell_val[i] = edge.back();
            ghosts[i] = (std::vector<unsigned int>{edge[0],edge[1]});
            i++;     
       }
    }
}

void mesh::sort_mesh()
{
    //std::cout << "Sorting mesh...\n";
    std::vector<std::vector<std::vector<int>>> mask =  {{{0,1}},
                                                        {{0,1},{1,0}},
                                                        {{0,1},{1,2},{2,0}},
                                                        {{0,1},{1,2},{2,3},{3,0}}}; // quad,trig
    std::vector<int> node_vec;

    vec2ui polygons = quads;

    polygons.insert(polygons.end(),trigs.begin(),trigs.end());
    polygons.insert(polygons.end(),ghosts.begin(),ghosts.end());

    int neighbour;

    int mask_i;
    int c_idx = 0;
    unsigned int wall_idx = 0;
    for(auto const& cell : cells)
    {
        mask_i = cell.N_faces-1;
        
        for(int w = 0; w < cell.N_faces; w++)
        {
            node_vec = {int(polygons[c_idx][mask[mask_i][w][0]]),
                        int(polygons[c_idx][mask[mask_i][w][1]])};


            neighbour = find_neigbour_cell(polygons, node_vec, c_idx);

            face wall(nodes[node_vec[0]],nodes[node_vec[1]], c_idx, neighbour,node_vec[0],node_vec[1]);
            wall.owner_cell_index = c_idx;
            wall.neigbour_cell_index = neighbour;

            if(wall_uniqueness(wall) && neighbour != -1)
            {
                walls.push_back(wall);   
                cells[c_idx].add_cell_wall(wall_idx);
                
                cells[neighbour].add_cell_wall(wall_idx);
                wall_idx++;
            }

            //std::cout << "cell n: " << c_idx << " has neighbour n: " << neighbour << "\n";
        }
        c_idx++;
    }

    //std::cout << "Mesh sorted...\n\n";
}

void mesh::set_owner_idx()
{
    int i = 0;
    for(auto& cell : cells)
    {
        int j = 0;
        //for(auto const& w : cell.cell_walls)
        for(int w = 0; w < cell.N_faces; w++)
        {
            if(walls[cell.cell_walls[w]].owner_cell_index == i)
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
    //std::cout << "Constructing cells...\n";

    std::vector<unsigned int> ghost_vec(3,0);

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
        cells[k+N_quads+N_trigs] = cell(1,ghosts[k],nodes);
    }

    //std::cout << "Done!\n\n";
}

void mesh::expand_ghost_cells()
{
    double x_g,y_g;
    std::vector<double> G = {0,0,0};
    ghost_node_offset = nodes.size();

    for(int i = N_cells; i < N; i++)
    {
        cells[walls[cells[i].cell_walls[0]].owner_cell_index];

        x_g = walls[cells[i].cell_walls[0]].xf + walls[cells[i].cell_walls[0]].n[0] * walls[cells[i].cell_walls[0]].S;
        y_g = walls[cells[i].cell_walls[0]].yf + walls[cells[i].cell_walls[0]].n[1] * walls[cells[i].cell_walls[0]].S;

        G = {x_g,y_g,0.0};

        nodes.push_back(G);

        ghosts[i-N_cells].push_back(nodes.size()-1);

        cells[i].ghost(ghosts[i-N_cells],nodes);
    }
}

void mesh::vertex_init()
{
    unsigned int idx;

    vec2ui polygons = quads;

    polygons.insert(polygons.end(),trigs.begin(),trigs.end());
    polygons.insert(polygons.end(),ghosts.begin(),ghosts.end());

    for(auto& wall : walls)
    {
        for(auto& vertex : wall.vertices)
        {
            idx = vertex.node_idx;
            //loop over all quads/trigs and find those containing certain node idx
            unsigned int cell_idx = 0;
            for(auto const& poly : polygons)
            {
                if(std::find(poly.begin(),poly.end(),idx) != poly.end())
                {
                    vertex.common_cell_idx.push_back(cell_idx);
                }
                cell_idx++;
            }

            vertex.common_cell_idx.shrink_to_fit();

            //check if found node is already stored in vertex struct
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

        if(!std::count(group_idx.begin(),group_idx.end(),val))
        {
            group_idx.push_back(val);
            group_size.push_back(0);
            boundary_groups.push_back(group());
        }

        
        for(uint i = 0; i < group_idx.size(); i++)
        {
            if(val == group_idx[i])
            {
                group_size[i]++;
                boundary_groups[i].member_idx.push_back(c + N_cells);
                boundary_groups[i].group_value = val;
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
}

void mesh::group_surfaces(std::vector<uint> data)
{
    uint val;
    std::vector<int> group_idx;
    std::vector<int> group_size;

    for(int c = 0; c < data.size(); c++)
    {
        val = data[c];

        if(!std::count(group_idx.begin(),group_idx.end(),val))
        {
            group_idx.push_back(val);
            group_size.push_back(0);
            physical_surface.push_back(group());
        }

        for(uint i = 0; i < group_idx.size(); i++)
        {
            if(val == group_idx[i])
            {
                group_size[i]++;
                physical_surface[i].member_idx.push_back(c);
                physical_surface[i].group_value = val;
            }
        }
    }

    int g = 0;
    std::cout << "Surface groups: \n";
    for(auto const& group : physical_surface)
    {
        std::cout << "group: " << g << " ,value: " << group.group_value
                  << " ,size: " << group.member_idx.size() << "\n";
        g++;
    }
}

void mesh::export_mesh()
{
    std::ofstream f(name + "_walls.txt");
    f << N_walls << "\n";

    int w_idx = 0;
    for(auto const& wall : walls)
    {
        f << w_idx << "\n";
        f << wall.xf << " " << wall.yf << " " << wall.S << "\n";
        f << wall.n[0] << " " << wall.n[1] << "\n";
        //f << wall.s[0] << " " << wall.s[1] << "\n";
        f << wall.owner_cell_index << " " << wall.neigbour_cell_index << "\n";
        f << "\n";
        w_idx++;
    }
    f.close();

    std::ofstream ff(name + "_cells.txt");
    ff << N_cells << " " << N_ghosts << "\n";
    int c_idx = 0;
    for(c_idx = 0; c_idx < N; c_idx++)
    {
        ff << c_idx << "\n";
        ff << cells[c_idx].x << " " << cells[c_idx].y << " " << cells[c_idx].V << "\n";

        ff << cells[c_idx].owner_idx[0] << " "
           << cells[c_idx].owner_idx[1] << " "
           << cells[c_idx].owner_idx[2] << " "
           << cells[c_idx].owner_idx[3] << "\n";

        for(auto const& wall : cells[c_idx].cell_walls)
        {
            ff << wall << " ";
        }
        ff << "\n";

        ff << "\n\n";
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
        // C.s[0] = double_vec[0];
        // C.s[1] = double_vec[1];

        std::getline(f, text);
        double_vec = extract(text);
        C.owner_cell_index = double_vec[0];
        C.neigbour_cell_index = double_vec[1];

        walls.push_back(C);
    }

    file.close();
}

void mesh::transform_axi_X()
{
    for(auto& cell : cells)
    {
        cell.V = cell.x*cell.V;
    }

    for(auto& wall : walls)
    {
        wall.S = wall.xf*wall.S;
    }
}

void mesh::transform_axi_Y()
{
    for(auto& cell : cells)
    {
        cell.V = (cell.y)*cell.V;
    }

    for(auto& wall : walls)
    {
        wall.S = (wall.yf)*wall.S;
    }
}

void mesh::transform_axisymetric(int axis_idx)
{
    std::cout << "Axisymetric transformation along ";

    switch (axis_idx)
    {
    case 0:
        std::cout << "x axis\n";
        transform_axi_X();
        break;
    case 1:
        std::cout << "y axis\n";
        transform_axi_Y();
        break;
    default:
        std::cout << "\n Error, axis index not available\nCombustor exiting...\n";
        exit(1);
    }


}