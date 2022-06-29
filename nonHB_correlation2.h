#ifndef NONHB_CORRELATION2_H
#define NONHB_CORRELATION2_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class nonHB_correlation2
{
    public:

    nonHB_correlation2();
    ~nonHB_correlation2();

    void Routine();
    void calc(CellFile &cel);
    void output();
    double** HB_start_time;
    //double** HB_corr;
    int* HB_nrecord;
    int** HB_water_index;
    int ndt;
    int ndt1;
    int nHB;
    ofstream ofs_nHB;
};

#endif