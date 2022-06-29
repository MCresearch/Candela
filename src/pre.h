#ifndef PRE_H
#define PRE_H

#include "cell.h"
#include "water.h"
#include "mj.h"


// Presolvation structure of ions 
class PRE
{
	public: 
	
	PRE();
	~PRE();

	void Routine();

	private:

	int count_geometry_number;
	int* ion_pt;
	int nion;
	ProtonTransfer PT;

	void search_presolvation(const Cell &cel, const int &igeo, int* ion_nacc, int* ion_ndon);
};


#endif
