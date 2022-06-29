#ifndef MOVIE_H
#define MOVIE_H

#include "cellFile.h"

class Movie 
{
public:
    
	Movie();
    ~Movie();
    
	void Routine();

	private:

	void snapshot(ofstream &ofs, const Cell &cel, const int &iat);
	void print(const Cell &cel, const int &it, const int &ia, ofstream &ofs, string &id);

	void print_special_hexane_only_ion(const Cell &cel, Water *water, ofstream &ofs, 
		const int &ito, const int &itc, const int &ith, const int &igeo);

	void print_special_hexane_all_water(const Cell &cel, Water *water, ofstream &ofs, 
		const int &itc, const int &igeo);
	void print_within_distance(const Cell &cel, const int &ia, ofstream &ofs, Water* &water, const int &igeo, const int &ito, const int &ith);
	double special_distance(const Vector3<double> &pos1, Vector3<double> &pos2, const double &celldm1, const double &celldm2, const double &celldm3);
};

#endif //Movie
