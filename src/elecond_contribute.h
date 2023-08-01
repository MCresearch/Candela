#ifndef ELECOND_CONTRIBUTE_H
#define ELECOND_CONTRIBUTE_H

#include "wfRead.h"
#include "binfstream.h"
#include "wavefunc.h"

#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma

/**
 * @brief calculate contribution of each band to conductivities
 * 
 */
class Elecond_contribute
{
  public: 
	
	Elecond_contribute(){};
	~Elecond_contribute(){};

	void Routine();

  private:
    void cal_contribute();
	

};


#endif