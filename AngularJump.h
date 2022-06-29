#ifndef ANGULARJUMP_H
#define ANGULARJUMP_H

#include "input.h"
#include "HBs.h"


class AngularJump
{
    public:
    static void setup_AJ(int nAJ, AngularJump* AJ);
    int AJ_ss_index;
    double AJ_ss_time;
    int centralO;
    int* Oab;
    int Oa;
    int Ob;
    double center_angle;
    int H;
    AngularJump();
    ~AngularJump();
};

#endif