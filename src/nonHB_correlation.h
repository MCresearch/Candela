#ifndef NONHB_CORRELATION_H
#define NONHB_CORRELATION_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class nonHB_correlation
{
    public:

    nonHB_correlation();
    ~nonHB_correlation();

    void Routine();
    void calc(CellFile &cel);
    void output();
    double*** nonHB_time;
    //double*** HB_corr;
    bool* inshell;
    bool* previous_inshell;
    int* nslice; 
    //int** HB_nrecord;
    int** HB_water_index;
    int ndt;
    int nHB;

    static bool HB_bonded(Water* water, int &iwater1, int &iwater2);
};

#endif
