#ifndef DOUBLE_DONOR_H
#define DOUBLE_DONOR_H

#include "input.h"
#include "HBs.h"


class DoubleDonor
{
    public:
    static void setup_DD(int &nDD, DoubleDonor* DD);
    double DD_start_time;
    double DD_end_time;
    int centralO;
    int* Oab;
    int Obefore;
    int Oafter;
    int H;
    string jump;
    DoubleDonor();
    ~DoubleDonor();
};

#endif