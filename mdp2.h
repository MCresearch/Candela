#ifndef MDP2_H
#define MDP2_H

#include "cell.h"

// mean density profile based on
// mean liquid interfaces
class MDP2
{
	public: 
	
	MDP2(){};
	~MDP2(){};

	void Routine();
	
	private:

	double* mean_density; // for average water distribution
	double dr;
	double rcut;
	int nmesh;

	double** interface;
	double*** gradient;

	void proximity(ifstream &ifs, const Cell &cel, const int &func);
	void read_write_ili(ifstream &ifs);

};


#endif
