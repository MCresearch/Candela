#ifndef EXTEND_H
#define EXTEND_H

#include "cell.h"

// each file containing a set of type of atoms.
class Extend
{
	public: 
	
	Extend(){};
	~Extend(){};

	static void Routine();
	
	private:

	static void ExtendCell( const Cell &cel_in, Cell &cel_out);

};


#endif
