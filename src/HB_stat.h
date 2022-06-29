#ifndef HB_STAT_H
#define HB_STAT_H

#include "HBs.h"
#include "cell.h"

class HB_stat
{
    public:
    HB_stat();
    ~HB_stat();

    void Routine();

    private:

    int* nHBs;
    int* HB_lifespan_hist;
    int* relative_vel_hist;
    int* angle_hist;
    int* r_oh_hist;
    int* r_oo_hist;

    double** HB_lifespan;
    double** incoming_time;
    //for func_b = 2
    double** last_bonded_time;

    double** angle;
    double** r_oh;
    double** r_oo;

    int** angle_r_oo_hist;
    int** angle_r_oh_hist;
    bool** accepted;
    void calc(const Cell &cel, const int &count_geometry_number,\
    int* &nHBs, int* &HB_lifespan_hist, int* &relative_vel_hist, int* &r_oo_hist, int* &r_oh_hist, int* &angle_hist, \
    double** &HB_lifespan, double** &incoming_time, double** last_bonded_time, bool** &accepted, double** &r_oo, double** &r_oh, double** &angle, 
    int** &angle_r_oo_hist, int** &angle_r_oh_hist);

};

#endif