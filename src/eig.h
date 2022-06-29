#ifndef EIG_H
#define EIG_H

#include "cell.h"

class Eig
{
public:
    
	Eig();
    ~Eig();
    
	//void Init();
	void Routine();

	private:

	void read_eig(Cell &cel, double* dis_eig); // distribution of eigenvalues.
	void snapshot(ifstream &ifs, ifstream &ifs2, const int &ic, Cell &cel);

	public:

	double* dis_eig;
	double dr;
	int nr;

};

#endif //eig
