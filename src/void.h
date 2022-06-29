#ifndef VOID_H
#define VOID_H

#include "cell.h"

class Void
{
	public: 
	
	Void(){};
	~Void(){};

	static void Routine();
	
	private:

	static void CreateVoids( const Cell &cel_in, Cell &cel_out);

};


#endif
