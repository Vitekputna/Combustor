#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "data_structures.h"
#include "initial_cond.h"
#include "boundary.h"


std::string read_between(std::string str, char from, char to)
{
    bool start = false, end = false;
    std::string res;

    for(auto const& c : str)
    {
        if(c == ';')
        {
            return "";
        }
        if(c == to)
        {
            end = true;
            return res;
        }

        if(start && !end)
        {
            res.push_back(c);
        }

        if(c == from)
        {
            start = true;
        }
    }

    return "";
}

std::vector<double> read_multiple_between(std::string str, char from, char del, char to)
{
    bool start = false, end = false;
    std::string word;
    std::vector<double> res;

    for(auto const& c : str)
    {
        if(c == ';')
        {
            return std::vector<double>{};
        }

        if(c == to)
        {
            //std::cout << word << "\n";
            res.push_back(std::stod(word));
            end = true;
            return res;
        }

        if(c == from)
        {
            start = true;
            continue;
        }

        if(c == del)
        {
            //std::cout << word << "\n";
            res.push_back(std::stod(word));
            word.clear();
            continue;
        }

        if(start && !end)
        {
            word.push_back(c);
        }
    }

    return std::vector<double>{};
}

std::string read_first_word(std::string str)
{
    std::string res;

    for(auto const& c : str)
    {
       if(c == ' ' || c == '=')
       {
            return res;
       }

       res.push_back(c);
    }

    return "";
}

int find_boundaryGroup(int bdr_val, boundary const& bdr, mesh const& msh)
{
    int i = 0;
    for(auto const& group : msh.boundary_groups)
    {
        if(bdr_val == group.group_value)
        {
            return i;
        }
        i++;
    }

    std::cout << "error: boundary with value: " << bdr_val << " not found\n";
    return -1;
}

void read_config_files(mesh const& msh, boundary& bdr, config& cfg, parameters& par, initial_conditions& IC)
{
    std::string str, keyword;
    double value;

    //initial data
    std::ifstream stream("config/init.txt");

    while(std::getline(stream,str))
    {
        if((str[0] != '/' && str[1] != '/') && !str.empty())
        {
            keyword = read_first_word(str);
            value = std::stod(read_between(str,'{','}'));

            if(keyword == "p_init")
            {
                IC.p_start = value;
            }
            else if(keyword == "T_init")
            {
                IC.T_start = value;
            }
            else if(keyword == "M_init")
            {
                IC.M_start = value;
            }
            else if(keyword == "alfa_init")
            {
                IC.alfa_start = value;
            }    
        }
    }

    //boundary data
    int g_idx;
    int bdrf_val;
    int mshGidx;
    std::vector<double> g_val;
    std::ifstream stream2("config/boundary.txt");

    while(std::getline(stream2,str))
    {
        if((str[0] != '/' && str[1] != '/') && !str.empty())
        {
            g_idx = std::stoi(read_between(str,'{','}'));
            g_val = read_multiple_between(str,'(',',',')');

            if(str.find("wall") != std::string::npos)
            {
                bdrf_val = 0;
                mshGidx = find_boundaryGroup(g_idx,bdr,msh);
                bdr.boundary_groups.push_back(boundary_group(bdrf_val,msh.boundary_groups[mshGidx].member_idx));

                //std::cout << g_idx << " " << bdrf_val - g_idx << "\n";
            }
            else if(str.find("supersonic_inlet") != std::string::npos)
            {
                bdrf_val = 1;
                mshGidx = find_boundaryGroup(g_idx,bdr,msh);
                bdr.boundary_groups.push_back(boundary_group(bdrf_val,msh.boundary_groups[mshGidx].member_idx));
                bdr.boundary_groups.back().p_0 = g_val[0];
                bdr.boundary_groups.back().T_0 = g_val[1];
                bdr.boundary_groups.back().Min = g_val[2];
                bdr.boundary_groups.back().alfa = g_val[3];

                supersonic_inlet(par,bdr.boundary_groups.back());

                //std::cout << g_idx << " " << bdrf_val - g_idx << "\n";
            }
            else if(str.find("supersonic_outlet") != std::string::npos)
            {
                bdrf_val = 2;
                mshGidx = find_boundaryGroup(g_idx,bdr,msh);
                bdr.boundary_groups.push_back(boundary_group(bdrf_val,msh.boundary_groups[mshGidx].member_idx));

                //std::cout << g_idx << " " << bdrf_val - g_idx << "\n";
            }
            else if(str.find("subsonic_inlet") != std::string::npos)
            {
                bdrf_val = 3;
                mshGidx = find_boundaryGroup(g_idx,bdr,msh);
                bdr.boundary_groups.push_back(boundary_group(bdrf_val,msh.boundary_groups[mshGidx].member_idx));                
                bdr.boundary_groups.back().p_0 = g_val[0];
                bdr.boundary_groups.back().T_0 = g_val[1];
                bdr.boundary_groups.back().alfa = g_val[2];

                subsonic_inlet(par,bdr.boundary_groups.back());

                //std::cout << g_idx << " " << bdrf_val - g_idx << "\n";
            }
            else if(str.find("subsonic_outlet") != std::string::npos)
            {
                bdrf_val = 4;
                mshGidx = find_boundaryGroup(g_idx,bdr,msh);
                bdr.boundary_groups.push_back(boundary_group(bdrf_val,msh.boundary_groups[mshGidx].member_idx));   
                bdr.boundary_groups.back().p_stat = g_val[0];

                subsonic_outlet(par,bdr.boundary_groups.back());
                
                //std::cout << g_idx << " " << bdrf_val - g_idx << "\n";
            }
        }
        
    }
    
    //controll data
    std::ifstream stream3("config/controll.txt");

    while(std::getline(stream3,str))
    {
        if((str[0] != '/' && str[1] != '/') && !str.empty())
        {
            keyword = read_first_word(str);
            value = std::stod(read_between(str,'{','}'));

            if(keyword == "min_iterations")
            {
                cfg.min_iter = value;
            }
            else if(keyword == "max_iterations")
            {
                cfg.max_iter = value;
            }
            else if(keyword == "max_residue")
            {
                cfg.max_res = value;
            }
            else if(keyword == "residue_idx")
            {
                cfg.res_idx = value;
            }
            else if(keyword == "bisection_iterations")
            {
                cfg.bisec_iter = value;
            }    
            else if(keyword == "timestep_freq")
            {
                cfg.n_t = value;
            }
            else if(keyword == "residue_freq")
            {
                cfg.n_r = value;
            }
            else if(keyword == "boundary_freq")
            {
                cfg.n_b = value;
            }
            else if(keyword == "export_freq")
            {
                cfg.n_exp = value;
            }
            else if(keyword == "CFL")
            {
                cfg.CFL = value;
            }
            else if(keyword == "timestep")
            {
                cfg.dt = value;
            }
        }
    }
}

