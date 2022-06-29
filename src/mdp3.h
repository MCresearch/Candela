#ifndef MDP3_H
#define MDP3_H

#include "cell.h"

// Copmute mean density profile (MDP) based on
// instantaneous liquid interfaces (ILI)
class MDP3
{
	public: 
	
	MDP3();
	~MDP3();

	void Routine();

	private:

	double** interface;
	double*** gradient;

	// two dimension data
	double** coord_xy;

	int nx;
	int ny;
	double x0;
	double y0;
	double dx;
	double dy;

	void proximity(ifstream &ifs, const Cell &cel);
};


#endif
