#ifndef OHO_ANGLE_H
#define OHO_ANGLE_H
#include "cellFile.h"
#include "HBs.h"

class oho_angle
{
public:
    oho_angle();
    ~oho_angle();

    void Routine();

    double* theta_distr;
    int ntheta;

    void calc(CellFile &cel, Water* water, int &ito, int &ith);
    void normalize();

};

#endif