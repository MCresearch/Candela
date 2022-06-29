#ifndef HB_STAT2_H
#define HB_STAT2_H

#include "HBs.h"
#include "cellFile.h"

class HB_stat2
{
    public:
    HB_stat2();
    ~HB_stat2();
    void Routine();
    void calc(const Cell &cel, int &nDD, int &nDA, int &nAA, double* &DD_distr, double* &DA_distr, double* &AA_distr);
    double* DD_distr_acc1;
    double* DD_distr_acc2;
    double* DD_distr_acc3;
    double* AA_distr_acc1;
    double* AA_distr_acc2;
    double* AA_distr_acc3;
    double* DA_distr_acc1;
    double* DA_distr_acc2;
    double* DA_distr_acc3;
    double* HH_distr;
    ofstream ofs_smallDD;
    ofstream ofs_smallDD_acc_don;
};
#endif