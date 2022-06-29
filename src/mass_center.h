#ifndef MASS_CENTER_H
#define MASS_CENTER_H

#include "atoms.h"

class mass_center
{
public:
	mass_center();
	~mass_center();

	void Routine();
	Vector3<double>* pre_wpos;
	Vector3<double>* wpos;
};
#endif