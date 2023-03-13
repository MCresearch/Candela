#ifndef Ele_Conductivity_H
#define Ele_Conductivity_H

#include "wfRead.h"
#include "binfstream.h"
#include "wavefunc.h"

#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma

//Calculate electric and thermal conductivity of electrons.
//By qianrui 2020-1-6
class Ele_Conductivity
{
	public: 
	
	Ele_Conductivity(){};
	~Ele_Conductivity(){};

	void Routine();

	void method1();
	void method2();

	private:
	void jjcorr_ks(const int ik, const int nt, const double dt,  WfRead& wfr,
            double* ct11, double* ct12, double* ct22);
	void calcondw(const int nt, const double dt, const double fwhmin, const double wcut, 
            const double dw_in, double *ct11, double *ct12, double *ct22);
	
	

};


#endif
/*
	calculation ele_conductivity
	wf_in_type   PWmat/QE1/QE2/ABACUS
	wfdirectory  path_to_wfs (if nscf==0, then we read wf in path_to_wfs; else if nscf=Nfiles(>0), then we read wf in path_to_multi_files/[1:Nfiles]/path_to_wfs)
	
	multi_directory	path_to_multi_files
	dw    0.01
	wcut  10
	smear  1
	fwhm 0.1
(or	n_fwhm 2
	fwhm 0.1 0.2	)
	temperature 1000 //K
	nscf   0/Nfiles
	
*/
