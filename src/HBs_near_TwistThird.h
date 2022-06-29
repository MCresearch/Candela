#ifndef HBS_NEAR_TWISTTHIRD_H
#define HBS_NEAR_TWISTTHIRD_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class HBs_near_TwistThird
{
    public:

    HBs_near_TwistThird();
    ~HBs_near_TwistThird();

    int** water_index;
    double** start_end_time;
    bool** HB_before_amid_after;

    void Routine();
    void allocate_index_time();
    void calc(CellFile &cel);
    bool Hbonded(Water* water, int &iwater1, int &iwater2);
};
#endif