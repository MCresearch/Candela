#ifndef TRAJADJ_H
#define TRAJADJ_H

#include "cell.h"
#include "water.h"


class Trajadj 
{
	public: 
	
	Trajadj();
	~Trajadj();

	void Routine();

	private:

	int itzone;
	void compute(const Cell &cel, double** trmap, const int &igeo);
	double dr; // delta r
	double rcut; // radius cutoff
	int nmesh; // number of mesh points.

};


#endif
