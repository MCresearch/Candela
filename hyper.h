#ifndef HYPER_H
#define HYPER_H

#include "cell.h"
#include "water.h"

// parameters for hyper coordinated 
class Hyper
{
	public: 
	
	Hyper();
	~Hyper();

	void Routine();
	static void dvec(const Cell &cel, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &da);

	private:

	void hypercoordinated(const Cell &cel, double* gr);

	double dr;
	double rcut;
	int nmesh;	
	int count_geometry_number;
	int igeo;

	double ave_donate;


	void four_accept(const Cell &cel, const Water* water, const int &ito, const int &ia, double* gr);
	void one_donate(const Cell &cel, const Water* water, const int &ito, const int &ia, double* gr);



};


#endif
