#ifndef HBS_CORRELATION_H
#define HBS_CORRELATION_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class HB_correlation
{
    public:

    HB_correlation();
    ~HB_correlation();

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