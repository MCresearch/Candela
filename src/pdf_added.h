#ifndef PDF_ADDED_H
#define PDF_ADDED_H

#include "cell.h"
#include "Honeycutt.h"

// mohan create this class on 2019-03-14
class PDF_ADDED
{
	public: 
	
	PDF_ADDED(){};
	~PDF_ADDED(){};

	static bool atom_in_ion(const Cell &cel, const int &it, const int &ia, int &io_of_ion);
	static bool atom_in_ion2(const Cell &cel, const int &it1, const int &ia, const int &it2, const int &ia2);
	static double compute_delta(const Cell &cel, const Water* water, const int &it_in, const int &ia);
	static bool water_ions(const Cell &cel, const Water* water, const int &io_of_ion, const int &ia2, bool &should_count);

};


#endif
