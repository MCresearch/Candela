#ifndef CationFreeEnergy_H
#define CationFreeEnergy_H
#include "cell.h"
#include "water.h"

class CationFreeEnergy
{
public:
    CationFreeEnergy(){};
    ~CationFreeEnergy(){};
    void Routine();
    double angle1(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3);
    double angle_vector(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3, Vector3<double> &vector1);
    double atom_vector(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &vector2);
    double angle3(double &v1x, double &v1y, double &v1z, double &v2x, double &v2y, double &v2z);
    double ** count_xy;

private:
    // convert pos to pdb
    void cal();
    int getAtomIndex(const Cell &cell, const int &ia1);
};

#endif