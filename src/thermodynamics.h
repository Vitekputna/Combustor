
inline double pressure(double gamma,double* U)
{
    return (gamma-1)*(U[3]-0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);
}