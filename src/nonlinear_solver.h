#pragma once
#include <iostream>
#include "boundary.h"

typedef unsigned int uint;

double bisection_method(double (*func)(boundary const&,double,double*),boundary const& B, double* params, int n_params) // params = [start,end,n_div,f_params......]
{
    double start = params[0];
    double end = params[1];
    double mid;
    unsigned int n_div = params[2];

    for(uint n = 0; n < n_div; n++)
    {
        mid = (start+end)/2;

        if(func(B,start,params+3)*func(B,mid,params+3) <= 0)
        {
            end = mid;
        }
        else
        {
            start = mid;
        }
    }

    //std::cout << mid << "\n"; 
    return (start+end)/2;
}