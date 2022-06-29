#ifndef ILI_H
#define ILI_H

#include "cell.h"

// instantaneous liquid interfaces
class ILI
{
	public: 
	
	ILI(){};
	~ILI(){};

	void Routine();
	
	private:

	double ref_rho;
	double** interface;
	double*** gradient;

	void get_ili(const Cell &cell);
	void get_ili2(const Cell &cell);
	double cal_gauss(const Cell &cel, const double &x_in, const double &y_in, const double &z_in, const bool &cal_grad,
	double* grad);

};


#endif
