#ifndef VELCOR_H
#define VELCOR_H

#include "cell.h"

//------------------------------------
// Velocity Correlation Function 
// Z(r) = <v1(t)v1(0)> / <v1(0)v1(0)>
//------------------------------------
class VelCor
{
	public: 
	
	VelCor(){};
	~VelCor(){};

	void Routine();
	
	private:

	void cal();

	float velocity_correlation_functions(
		const Cell &cel1,
		const Cell &cel2,
		const int &ia);

	void write_vc(
		float** vc,
		const int &nvel, 
		const bool* file_exit) const;

};


#endif
