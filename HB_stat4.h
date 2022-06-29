#ifndef HB_STAT4_H
#define HB_STAT4_H

#include "HBs.h"
#include "cellFile.h"

class HB_stat4
{
    public:
    HB_stat4();
    ~HB_stat4();
    void Routine();

    void calc(CellFile &cel);
    ofstream* ofs;
    ofstream ofs_DD;
    int nspace;
    double* DD_start_time;
    int* previous_donateO;
    int** previous_doubledonateO;
    int* now_donateO;
    bool* previous_double_donate;
    string* DD_HB;
    bool* double_donate;
    int ndouble_donate;
    int nHB_jump;
};

#endif