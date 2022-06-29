#ifndef ORIENTATION_TCF_H
#define ORIENTATION_TCF_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class Orientation_TCF
{
    public:

    Orientation_TCF();
    ~Orientation_TCF();

    void Routine();
    void calc(CellFile &cel, int &count_geometry_number);
    void output();
    void calc_tcf();
    //double** HB_start_time;
    //double** HB_corr;
    double* time_serial;
    Vector3<double>** orient_vec;
    int ndt;
    int ndt1;
    int n_mlc;
    //ofstream ofs_nHB;
};

#endif