#ifndef TUNE_STRU_H
#define TUNE_STRU_H

#include "cell.h"
#include "water.h"

// Tune structures
class Tune_Stru
{
	public: 
	
	Tune_Stru();
	~Tune_Stru();

	void Routine();

	private:

	void adjust_structure(const Cell &cell, const int &igeo, 
		ofstream &ofs, ofstream &ofs2);

	int count_geometry_number;

};


#endif
