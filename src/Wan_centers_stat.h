#ifndef WAN_CENTERS_STAT_H
#define WAN_CENTERS_STAT_H

#include "HBs.h"
#include "cellFile.h"

class Wan_centers_stat
{
    public:
    Wan_centers_stat();
    ~Wan_centers_stat();

    void Routine();
    static void allocate_wan(const Water* water, const Cell &cel, int** &bond_wan_index, int** &lone_wan_index, const int &ito, const int &ith);
    
    private:
    void Wan_angle_distr(Cell &cel);
    void calculate(const Water* water, const Cell &cel, int** &bond_wan_index, int** &lone_wan_index, const int &ito, const int &ith);
    void out(int &count_geometry_number);
    void put_back_cell(Vector3<double> &pos);

    int nangle;
    int nr;
    double **lone_2D_distr;
    double *lone_angle_distr;
};

#endif