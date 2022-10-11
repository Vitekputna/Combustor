#pragma once
#include <vector>
#include <string>
#include <limits>

typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<unsigned int>> vec2ui;
typedef std::vector<double> vec1d;
typedef std::vector<int> vec1i;
typedef std::vector<unsigned int> vec1ui;

typedef unsigned int uint;

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
    char N_faces;
    double x = 0, y = 0, V;
    std::vector<double> X,Y;
    std::vector<unsigned int> cell_walls;
    unsigned char free_wall_slot_idx = 0;
    int owner_idx[4] = {-1,-1,-1,-1};

    cell();
    cell(int N_walls, vec1ui nodes,vec2d const& all_nodes);
    void add_cell_wall(unsigned int wall_idx);
};

struct group
{
    group();
    uint group_value;
    std::vector<uint> member_idx;
};

typedef std::vector<cell> cell_vec;
typedef std::vector<face> face_vec;

class mesh
{
    public:
    std::string name;
    vec2d nodes;
    vec2ui edges,quads,trigs;

    cell_vec cells;
    face_vec walls;

    int quads_offset;
    vec2ui ghosts;
    vec1i ghost_cell_val;

    std::vector<group> boundary_groups;

    unsigned int N_cells, N_walls, N_ghosts, N_trigs, N_quads, N;
    double min_V = std::numeric_limits<double>::max();
    
    mesh(int i);
    mesh(std::string path);
    void load_mesh(std::string path, vec2d& nodes, vec2ui& edges, vec2ui& quads);

    template<typename T>
    std::vector<std::vector<T>> read_segment(std::vector<std::string>& text,int from,int offset, int length);

    void print_mesh();
    void sort_mesh();

    int find_neigbour_cell(vec2ui const& polygons, vec1i common_nodes_idx,int owner_idx);
    bool wall_uniqueness(face const& new_wall);
    void construct_ghost_cells();
    void construct_cells();
    void set_owner_idx();

    void group_inlets();
 
    std::vector<double> extract(std::string& text);

    void export_mesh();
    void import_mesh(std::string path);

    void transform_axisymetric(int axis_idx);

    private:

    void transform_axi_X();
    void transform_axi_Y();
};