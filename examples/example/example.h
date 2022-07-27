#ifndef EXAMPLE_H
#define EXAMPLE_H
#include "Cellfile.h"
class Example
{
    public:
    Example();
    ~Example();

    void Routine();
    void calc(CellFile &cel);
    void output();

    double* example_arr;
    int example_len;
};

#endif