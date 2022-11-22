#pragma once
#include <vector>
#include "data_structures.h"

struct initial_conditions
{
    initial_conditions(int n_comp)
    {
        composition_mass_frac.resize(n_comp,0.0);
    }
    std::vector<double> U;
    // std::vector<int> composition;
    std::vector<double> composition_mass_frac;
    double p_0, T_0, Min, alfa, beta, p_stat, p_start, M_start, T_start, alfa_start;
};