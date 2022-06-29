#ifndef HBS_NEAR_DOUBLEDONOR_H
#define HBS_NEAR_DOUBLEDONOR_H

#include "input.h"
#include "HBs.h"
#include "DoubleDonor.h"
#include "cellFile.h"
#include "dist2.h"

class HBs_near_DoubleDonor
{
    public:

    HBs_near_DoubleDonor();
    ~HBs_near_DoubleDonor();

    void Routine();
    //void double_array_new(double** arr, int dim1, int dim2);
    double** OaOH_angle;
    double** ObOH_angle;
    double** OOO_angle;
    double** OH_central_angle;
    double** OH_plane_angle;

    double** OaO_distance;
    double** ObO_distance;
    
    double** central_anthH_donate;
    double** central_accept;
    double** donOa_donate;
    double** donOb_donate;
    double** donOa_accept;
    double** donOb_accept;
    DoubleDonor* DD;

    void calc(CellFile &cel);
    int ndt;
    void output();
    void single_output(double** arr, string jump_file, string rattle_file);
};

#endif