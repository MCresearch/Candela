#ifndef IN_PDF3_H
#define IN_PDF3_H

#include "cell.h"
#include "pdf_added.h"
#include "water.h"
class incrementalPDF3
{
public:
	
	incrementalPDF3();
	~incrementalPDF3();

	void Routine();

private:
    int* one_path;
    // 0: DDD/AAA; 1: DDA/DAA; 2: DAD/ADA; 3: ADD/AAD
    int* two_paths;
    // 0: 00; 1: 01; 2: 02; 3: 03; 4: 11; 5: 12; 6: 13; 7: 22; 8: 23; 9: 33
    int* npaths;
    bool* calculated;

    double** gr_one_path;
    double** gr_two_paths;
    double* gr_three_paths;

    void calc(Cell &cel);
    void allocate();
    void set_false(Cell &cel, Water *water, int &iwater);

    int* iwater_in_shell;
    int nwater_in_shell;

    int* type_of_path;
    int npaths_to_water;

    void search_third_shell_water(Cell &cel, Water *water, int &iwater, int &ito, int &ith);
    void set_zero();
    void search_paths(Cell &cel, Water *water, int &iwater, int &iwater_select, int &ito, int &ith);
    void out(int &count_geometry_number);
};

#endif