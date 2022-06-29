#ifndef MATH_H
#define MATH_H

#include "gfun.h"

class Math
{
public:

    Math(){};
    ~Math(){};

    static void Simpson_Integral
    (
        const int mesh,
        const double *func,
        const double *rab,
        double &asum
    );
};

#endif
