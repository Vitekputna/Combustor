#pragma once

inline void no_source_cartesian(double dt,double* W){}

inline void heat_source_cartesian(double dt,double* W)
{
    W[3] += dt*3e6;
}