#pragma once
#include <vector>
#include <string>
#include <limits>

typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<unsigned int>> vec2ui;
typedef std::vector<double> vec1d;
typedef std::vector<int> vec1i;
typedef std::vector<unsigned int> vec1ui;


class face
{
    public:
    double xf,yf,S;
    double n[2], s[2];
    int owner_cell_index, neigbour_cell_index;
    face(vec1d const& a, vec1d const& b);
    face();
};

class cell
{
    public:
    double x = 0, y = 0, V;
    unsigned int cell_walls[4];
    unsigned char free_wall_slot_idx = 0;
    int owner_idx[4] = {-1,-1,-1,-1};

    cell();
    cell(vec1ui nodes,vec2d const& all_nodes);
    void add_cell_wall(unsigned int wall_idx);
};

typedef std::vector<cell> cell_vec;
typedef std::vector<face> face_vec;

class mesh
{
    public:
    std::string name;
    vec2d nodes;
    vec2ui edges,quads;
    cell_vec cells;
    face_vec walls;

    int quads_offset;
    vec1i ghost_cell_idx;
    vec1i ghost_cell_val;

    unsigned int N_cells, N_walls, N_ghosts, N;
    double min_V = std::numeric_limits<double>::max();
    
    mesh(int i);
    mesh(std::string path);
    void load_mesh(std::string path, vec2d& nodes, vec2ui& edges, vec2ui& quads);
    void print_mesh();
    void sort_mesh();

    int find_neigbour_cell(vec1i common_nodes_idx,int owner_idx);
    bool wall_uniqueness(face const& new_wall);
    void construct_ghost_cells();
    void construct_cells();
    void set_owner_idx();

    void group_inlets();
 
    std::vector<double> extract(std::string& text);

    void export_mesh();
    void import_mesh(std::string path);

    
};