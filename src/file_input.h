#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include "data_structures.h"
#include "source_functions.h"
#include "initial_func.h"
#include "boundary.h"
#include "solver.h"

// predělat tohle na volání konfiguračních funkcí a ne if elseif špagetu


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

std::vector<std::vector<double>> read_config_files(mesh& msh, boundary& bdr, config& cfg, parameters& par, solver& sol)
{
    std::string str, keyword;
    double value;

    //controll data
    std::ifstream stream3("config/controll.txt");

    while(std::getline(stream3,str))
    {
        if((str[0] != '/' && str[1] != '/') && !str.empty())
        {
            keyword = read_first_word(str);
            

            if(keyword == "min_iterations")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.min_iter = value;
            }
            else if(keyword == "max_iterations")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.max_iter = value;
            }
            else if(keyword == "max_residue")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.max_res = value;
            }
            else if(keyword == "residue_idx")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.res_idx = value;
            }
            else if(keyword == "bisection_iterations")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.bisec_iter = value;
            }    
            else if(keyword == "timestep_freq")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.n_t = value;
            }
            else if(keyword == "residue_freq")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.n_r = value;
            }
            else if(keyword == "boundary_freq")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.n_b = value;
            }
            else if(keyword == "export_freq")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.n_exp = value;
            }
            else if(keyword == "CFL")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.CFL = value;
            }
            else if(keyword == "timestep")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.dt = value;
            }
            else if(keyword == "max_time")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.max_time = value;
            }
            else if(keyword == "export_interval")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.export_interval = value;
            }
            else if(keyword == "axisymetric")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.r_variable_idx = (int)(value);   
                
                sol.source_func = axisymetric_source;
                sol.flux_func = HLL_flux_axi;

                msh.transform_axisymetric(value);
            }
            else if(keyword == "dimension")
            {
                value = std::stod(read_between(str,'{','}'));
                cfg.dim = (int)(value+2);
                cfg.vel_comp = value;
            }


        }
    }

    //initial data
    std::ifstream stream("config/init.txt");

    int sur_idx = 0;
    int n_ps = 0;

    std::vector<initial_conditions> IC_vec;

    while(std::getline(stream,str))
    {
        if((str[0] != '/' && str[1] != '/') && !str.empty())
        {
            keyword = read_first_word(str);
            value = std::stod(read_between(str,'{','}'));

            if(keyword == "surface")
            {
                sur_idx = value;
                IC_vec.resize(sur_idx+1);
                n_ps++;
            }

            if(keyword == "p_init")
            {
                IC_vec[sur_idx].p_start = value;
            }
            else if(keyword == "T_init")
            {
                IC_vec[sur_idx].T_start = value;
            }
            else if(keyword == "M_init")
            {
                IC_vec[sur_idx].M_start = value;
            }
            else if(keyword == "alfa_init")
            {
                IC_vec[sur_idx].alfa_start = value;
            }    
            else if(keyword == "beta_init")
            {
                IC_vec[sur_idx].beta = value;
            }    
        }
    }

    if(n_ps != IC_vec.size())
    {
        std::cout << "Error creating physical surface BC, check your settings...\n";
        exit(1);
    }

    std::vector<std::vector<double>> ret_vec;

    for(int i = 0; i < IC_vec.size(); i++)
    {
        ret_vec.push_back(move_flow(cfg.dim,IC_vec[i]));
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
                bdr.boundary_groups.back().beta = g_val[4];

                supersonic_inlet(cfg.dim,par,bdr.boundary_groups.back());

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
                bdr.boundary_groups.back().beta = g_val[3];

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

    return ret_vec;
}


