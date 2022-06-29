#ifndef VACUUM_H
#define VACUUM_H

#include "cell.h"

// each file containing a set of type of atoms.
class Vacuum
{
	public: 
	
	Vacuum(){};
	~Vacuum(){};

	static void Routine();
	
	private:

	static void AddVacuum( const Cell &cel_in, Cell &cel_out);

};


#endif
