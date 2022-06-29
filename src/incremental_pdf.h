#ifndef IN_PDF_H
#define IN_PDF_H

#include "cell.h"
#include "pdf_added.h"

class incrementalPDF
{
public:
	
	incrementalPDF();
	~incrementalPDF();

	void Routine();
	void bubble_sort(double* list, const int &length);
	void bubble_sort2(double* list, int* index_list, const int &length);
	void topological_neighbor_PDF(CellFile &cel, Water* &water, double** &gr_tmp, int* &total_nwater_in_shell, int nshell, int ito, int ith);

	double dr;
	double rcut;
	int nmesh;
	double** gr;
	double** gr_tmp;
	double** distance_matrix;
	int** index_matrix;
	int ntheta;
	int nq;
	double* tetraq_distr;
	double*** theta_r_oo;
	double*** tetraq_r_oo;
	double rho_ion;
};

#endif
