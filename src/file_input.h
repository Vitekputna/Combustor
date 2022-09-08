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

void read_config_files(boundary& bdr, config& cfg, parameters& par, initial_conditions& IC)
{
    std::string str, keyword;
    double value;

    //initial data
    std::ifstream stream("config/init.txt");
    while(std::getline(stream,str))
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
        //std::cout << keyword << " " << value << "\n";
    }

    //boundary data
    int g_idx;
    int bdrf_val;
    std::vector<double> g_val;
    std::ifstream stream2("config/boundary.txt");
    while(std::getline(stream2,str))
    {
        g_idx = std::stoi(read_between(str,'{','}'));
        g_val = read_multiple_between(str,'(',',',')');

        if(str.find("wall") != std::string::npos)
        {
            bdrf_val = 0;
            bdr.boundary_func_mask.push_back(bdrf_val - g_idx);
        }
        else if(str.find("supersonic_inlet") != std::string::npos)
        {
            IC.p_0 = g_val[0];
            IC.T_0 = g_val[1];
            IC.Min = g_val[2];
            IC.alfa = g_val[3];

            bdrf_val = 1;

            bdr.boundary_func_mask.push_back(bdrf_val - g_idx);
        }
        else if(str.find("supersonic_outlet") != std::string::npos)
        {
            bdrf_val = 2;

            bdr.boundary_func_mask.push_back(bdrf_val - g_idx);
        }
        else if(str.find("subsonic_inlet") != std::string::npos)
        {
            IC.p_0 = g_val[0];
            IC.T_0 = g_val[1];
            IC.alfa = g_val[2];

            bdrf_val = 3;

            bdr.boundary_func_mask.push_back(bdrf_val - g_idx);
        }
        else if(str.find("subsonic_outlet") != std::string::npos)
        {
            IC.p_stat = g_val[0];

            bdrf_val = 4;

            bdr.boundary_func_mask.push_back(bdrf_val - g_idx);
        }
        
    }
    


}

