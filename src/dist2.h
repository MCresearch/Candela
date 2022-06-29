#ifndef DIST2_H
#define DIST2_H

#include "cell.h"
#include "water.h"

// 3D distribution
class Dist2
{
	public: 
	
	Dist2();
	~Dist2();
    static void putback_cell(const Vector3<double> &pos1, Vector3<double> &pos2);

	void Routine();

    private:
    double dr;
    double*** dist3D;
    double*** dist3D2;
    int count_geometry_number;

    double* which_HOO_angle;

    void calculate_3D(Cell &cel);
    bool setup_axis(const int &iwater, const int &iH, const int &ito, const int &ith, const Cell &cel, const Water &water, Vector3<double> &xaxis, Vector3<double> &yaxis, Vector3<double> &zaxis);
    void wannier_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, const Water &water,
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, int** bond_wan_index, int** lone_wan_index);

    int nabnormal;
    bool fill_shell(const Cell &cel, const Water* water, const int &iwater, 
    int* &first_shell_index, int** &second_shell_index, int &nwater_first_shell, int* &nwater_second_shell);

    void second_shell_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water,
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH, 
    int* &first_shell_index, int** &second_shell_index, int &nwater_first_shell, int* &nwater_second_shell);
    void first_shell_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water, 
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH);
    void second_shell_3D2(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water, 
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH);
    void calc_which_HOO_angle(const Cell &cel, Water* &water, const int &ito, const int &ith, const int &iH, const int &water_index, int &water_index1, int &water_index2);
    void addup_3D(const Cell &cel, const int &ito, const int &iwater, const int &Oindex, const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis);
};
#endif
