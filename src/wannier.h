#ifndef WANNIER_H
#define WANNIER_H

#include "cell.h"

class Wannier
{
public:
    
	Wannier();
    ~Wannier();
    
	void Routine(); // compute dipole moments
	// output "all_dipole.dat"
	// output "all_vdipole.dat" 
	// output "total_dipoleM.dat"
	// output "distribution_dipole.dat"
	// output "distribution_MLWF.dat"


	void infrared(); // compute infrared spectra based on dipole moments

	private:

	void correlation_function_vdipole(double* uacf); // calculate the correlation function of vdipole
	void correlation_function_dipole(double* uacf); // calculate the correlation funtion of dipole
	void read_correlation_vdipole(double* uacf); // read the correlation function of vdipole
	void read_correlation_dipole(double* uacf); // read the correlation function of dipole
	void fourier_transform(double* uacf, const int &ndim); // perform the Fourier transform
	int tcor; // number of snapshots for correlation function (correlation length)
	double interval;
	double dt_au; // dt in atomic unit

	private:

	void read_wan(Cell &cel, double* dis_Hwan, ofstream &ofs_dipole, ofstream &ofs_vdipole, ofstream &ofs_tdipole);
	void snapshot(ifstream &ifs, ifstream &ifs2, const int &ic, Cell &cel);

	public:

	double* dis_Hwan;
	double dr;
	int nr;
	double** coord_xy;

	// store previous dipole in order to compute
	// the velocities of dipole moments
	bool init_previous;
	double previous_t;
	double* previous_d;
	double* previous_dx; // t(i-1)
	double* previous_dy;
	double* previous_dz;
	double previous2_t;
	double* previous2_d; // t(i-2)
	double* previous2_dx;
	double* previous2_dy;
	double* previous2_dz;

	// dipole
	int nd;
	double* dipole_dis;
	double avg_dipole;
	int count_dipole;

	// total dipole
	double rcut_total_dipole;
	double dr_total_dipole;
	double* dis_total_dipole;

};

#endif //Wannier
