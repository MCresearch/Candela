#ifndef INSERT_H
#define INSERT_H

#include "cell.h"

// each file containing a set of type of atoms.
class Insert
{
	public: 
	
	Insert(){};
	~Insert(){};

	static void Routine();
	
	private:

	static void InsertAtoms( const Cell &cel_in, Cell &cel_out);

};


#endif
