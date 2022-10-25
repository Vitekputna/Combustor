#pragma once
#include <vector>
#include "data_structures.h"

struct initial_conditions
{
    std::vector<double> U;
    double p_0, T_0, Min, alfa, beta, p_stat, p_start, M_start, T_start, alfa_start;
    parameters par;
};