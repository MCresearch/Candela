#ifndef VEL_H
#define VEL_H

#include "cell.h"

class Vel
{
	public: 
	
	Vel(){};
	~Vel(){};

	void Routine();
	
	private:

	void cal( Cell &cel_in );


};


#endif
