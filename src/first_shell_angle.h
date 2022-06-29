#ifndef FIRST_SHELL_ANGLE_H
#define FIRST_SHELL_ANGLE_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class First_shell_angle
{
    public:

    First_shell_angle();
    ~First_shell_angle();

    void Routine();
    void calc(CellFile &cel);
    void output();

    double* D_HOO_angle;
    double* A_OOc_angle;
    double* A_alpha_angle;
    double* A_beta_angle;

    double* D_alpha_angle;
    double* D_beta_angle;

    int nangle;
};

#endif
