#include "data_structures.h"

double time_step(variables& var)
{
    for(int i = 0; i < var.N; i++)
    {
        var.W(i,0);
    }

    return 0;
}
