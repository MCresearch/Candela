#ifndef HB_STAT3_H
#define HB_STAT3_H

#include "HBs.h"
#include "cellFile.h"

class HB_stat3
{
    public:
    HB_stat3();
    ~HB_stat3();
    void Routine();
    double* OOH_angle;
    double* OOO_angle;
    double* OH_central_angle;
    double* OH_plane_angle;

    double* OO_distance;
    
    double* central_donate;
    double* central_accept;
    double* donO_donate;
    double* donO_accept;
    double** don_acc;
    double ndouble_donate;
    ofstream ofs_mid;

    void calc(Cell &cel);
    double calc_OH_central_plane_angle(Cell &cel, Water* water, const int &iwater, const int &iwater1, const int &iwater2, const int &iH, const int &ito, const int &ith);
    void normalize_output(int &count_geometry_number);
    static void single_normalize(double* arr, double &rcut, double &dr);
    static void single_output(double* arr, double &rcut, double &dr, string file_name, bool half_plus);
};
#endif