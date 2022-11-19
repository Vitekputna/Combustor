#pragma once
#include <vector>
#include "data_structures.h"

struct initial_conditions
{
    std::vector<double> U;
    std::vector<int> composition;
    std::vector<double> composition_mass_frac;
    double p_0, T_0, Min, alfa, beta, p_stat, p_start, M_start, T_start, alfa_start;
};