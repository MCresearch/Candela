#ifndef HONEYCUTT_H
#define HONEYCUTT_H

#include "cell.h"

//----------------------------------
// Honeycutt Anderson 
//----------------------------------
class Honeycutt
{
	public: 

	Honeycutt();
	~Honeycutt();

	void Routine();

	static void setup_nadj(const Cell &cel, int* nadj_in, int **adj_index_in, const int &max_num_adj_in);

	// nsn: number of shared neighbours
	// nsb: number of shared bonds
	static bool search_neighbours(const int &iat, 
				const int &iat2, int* nadj_in, int** adj_index_in, const int &max_num_adj_in,
				const int &nsn, const int &nsb);

	int max_num_adj;
	int* nadj;
	int** adj_index;

	private:

	void cal();
	void HA_pairs(const Cell &cel);
	void shared_neighbours(const int &iat, const int &iat2, int* nadj_in, int **adj_index_in, const int &max_num_adj_in);

	string *HA_name;
	double* HA_percent;
	int max_name;

	ofstream ofs_snapshot;
	bool* selected_atom;

};


#endif
