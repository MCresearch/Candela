#ifndef HBS_NEAR_ANGULARJUMP_H
#define HBS_NEAR_ANGULARJUMP_H

#include "input.h"
#include "HBs.h"
#include "AngularJump.h"
#include "cellFile.h"
#include "dist2.h"

class HBs_near_AngularJump
{
    public:

    HBs_near_AngularJump();
    ~HBs_near_AngularJump();

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
    AngularJump* AJ;

    void calc(CellFile &cel);
    int ndt;
    static void calc_OH_central_plane_angle(Cell &cel, Water* water, const int &iwater, const int &iwater1, const int &iwater2,  const int &iH, const int &ito, const int &ith, double &central_angle, double &planar_angle);
    void output();
    void single_output(double** arr, int dim1, int dim2, string file_name);
};

#endif