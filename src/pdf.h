#ifndef PDF_H
#define PDF_H

#include "cell.h"
#include "pdf_added.h"
#include "Honeycutt.h"
#include "input.h"
#include "HBs.h"
class PDF
{
	public: 
	
	PDF();
	~PDF();

	void Routine();
	
	private:

	// calculate the pair distribution function.
	void cal();

	// calculate the paris under periodic boudary conditions.
	void periodic_pairs( const Cell &cell, double * gr, const int &option, const int &index_pdf);
	void sort_multi_pdf(double** multi_pdf, int* npdf_count);
	// calculate the pairs for each (i+1)dr.
    void pairs( const Cell &cel, double *gr);
	void out_multiple_pdf(double** multi_pdf);
	bool correct_ion_correct_time(const Cell &cel, Water* water, int &ito, int io_of_ion);
	// useless
	//void static_structure_factor( const Cell &cel, const double &rho_ion, const double *gr, double *sf) const;
	
	private:
	ofstream ofs_calc_snapshots;
	double dr; // delta r
	double rcut; // radius cutoff
	int nmesh; // number of mesh points.
	int count_geometry_number;

	double rho_ion;
	double dg;
	int ng;
	void static_structure_factor(const double &rho_ion, const double *gr, double *sf) const;

	// for calculating std for 
	int npdf;
	double** multi_pdf;
	int* npdf_count;

	// for calculating pdf immediately in advance of PT
	double* snapshot_time_pt;
	int* iindex;
	int* iindex_p; 
	//static ifstream ifs_trans;
	bool* return_jump;


	Honeycutt HA; // mohan added 2017-01-15

};


#endif
