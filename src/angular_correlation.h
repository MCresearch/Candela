#ifndef ANGULAR_CORRELATION_H
#define ANGULAR_CORRELATION_H

#include "input.h"
#include "HBs.h"
#include "cellFile.h"
#include "dist2.h"

class angular_correlation
{
    public:

    angular_correlation();
    ~angular_correlation();

    void Routine();
    void calc(CellFile &cel);
    void output();
    Vector3<double>** OH_vector;
    int ngeo;
    int nrecord;
    //double** HB_corr;
};

#endif