#ifndef NONHB_CORRELATION3_H
#define NONHB_CORRELATION3_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class nonHB_correlation3
{
    public:

    nonHB_correlation3();
    ~nonHB_correlation3();

    void Routine();
    void calc(CellFile &cel);
    void output();
    double** nonHB_start_time;
    double** HB_start_time;
    //double** HB_corr;
    int* nonHB_nrecord;
    int* HB_nrecord;
    int** HB_water_index;
    int ndt;
    int ndt1;
    int nHB;
};

#endif