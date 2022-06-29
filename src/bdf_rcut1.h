#ifndef BDF_RCUT1_H
#define BDF_RCUT1_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class BDF_rcut1
{
    public:

    BDF_rcut1();
    ~BDF_rcut1();

    void Routine();
    void calc(CellFile &cel);
    void output();

    double* bdf_angle;

};

#endif
