#ifndef IN_PDF2_H
#define IN_PDF2_H

#include "cell.h"
#include "pdf_added.h"
#include "water.h"
class incrementalPDF2
{
public:
	
	incrementalPDF2();
	~incrementalPDF2();

	void Routine();
    static void put_back(Vector3<double> &pos1, Vector3<double> &pos2);
private:
    double** gr;
    void calc(Cell &cel);
    void allocate();
    void out();
    bool** calculated;
    void set_false(Cell &cel, Water* water, int &iwater);

    double calc_angle_planar(Cell &cel, Water* water, int &iwater, int &iwater2, int &iwater3, int &iwater4, int &ito, int &ith);
    
    double* angle_planar_total;
    double** angle_planar;
    double*** angle_planar_OOO;
    double** angle_planar_OOO_total;
    double*** angle_OOO_2; 
    double** angle_OOO_2_total;
    // This is to evaluate the joint distribution of angle O1-O2-O3 and angle O2-O3-O4.
    double*** distance_12_34_14;
    double** distance_12_34_14_total;

    ofstream ofs_xyz1;
    ofstream ofs_xyz2;
    ofstream ofs_xyz3;
    ofstream ofs_xyz4;

    int ngeometry;
    void output_geometry(Cell &cel, Water* water, int &iwater, int &iwater2, int &iwater3, int &iwater4, int &ito, int &ith, const int &type);
};

#endif