#ifndef HBS_CORRELATION2_H
#define HBS_CORRELATION2_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class HB_correlation2
{
    public:

    HB_correlation2();
    ~HB_correlation2();

    void Routine();
    void calc(CellFile &cel);
    void output();
    double*** HB_start_time;
    //double*** HB_corr;
    bool* HB_connected;
    bool* HB_previous_connected;
    int* nslice; 
    //int** HB_nrecord;
    int** HB_water_index;
    int ndt;
    int nHB;
    //ofstream ofs_nHB;
};

#endif