#ifndef REORGANIZE_H
#define REORGANIZE_H

#include "cell.h"

class Reorganize 
{
public:
    
	Reorganize();
    ~Reorganize();
    
	void Routine(); //reorganize the cell, geometry, and MLWF files 

	private:

	void print_out(Cell &cel, ofstream &ofsp, ofstream &ofsw, ofstream &ofsc, ofstream &ofsx);

};

#endif //Reorganize
