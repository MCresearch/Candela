#ifndef WANNIER1_H
#define WANNIER1_H

#include "cell.h"

class Wannier1
{
public:
    
	Wannier1();
    ~Wannier1();
    
	void Routine(); 

	private:

	void read_wan(Cell &cel, double* dis_Hwan);
	void angle_analysis(Cell &cel, double* dis_Hwan);
	double angle(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3);
	void analysis2(Cell &cel);
	void analysis3(Cell &cel);

	double dr;
	int nr;
	double* dis_Hwan;

	// two dimension data
	double** coord_xy;
	double*** data_xyz;
	
	int nx;
	int ny;
	int nz;
	double x0;
	double y0;
	double z0;
	double dx;
	double dy;
	double dz;

	public:
};

#endif //Wannier1
