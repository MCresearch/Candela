#ifndef BONDED_H
#define BONDED_H

#include "gfun.h"
#include "atom.h"

class Bonded
{
	public:
    
	Bonded(){};
    ~Bonded(){};

	void Init(Atom* atom, int* total_ring_number, int* ringtot_number, int &frame, double &time);

	private:

};

#endif //BONDED
