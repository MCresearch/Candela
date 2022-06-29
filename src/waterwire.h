#ifndef WATERWIRE_H
#define WATERWIRE_H

#include "cell.h"
#include "water.h"
#include "mj.h"


// Properties of Waterwire 
class Waterwire
{
	public: 
	
	Waterwire();
	~Waterwire();

	void Routine();

	private:

	int count_geometry_number;
	ProtonTransfer PT;

	void search_compression(const Cell &cel, const int &igeo);
	void single_jump(const Cell &cel, Water *water, const int &ito, const int &ith);
	void double_jump(const Cell &cel, Water *water, const int &ito, const int &ith);

	double** coord_xy;
	double* bond_df;
	double* bond_df1;
	double* bond_df2;
	double* dis_o1;
	double* dis_o2;
	double* dis_o3;
	double* dis_o123;
	int total_frames;
    int total_frames2;

	double avg_bonds; // bond angle ooo
	double avg_bonds1; // bond angle h1oo
	double avg_bonds2; // bond angle h2oo
	double avg_dis1;
	double avg_dis2;
	double avg_dis12;
	int count_bonds;
	int count_bonds1;
	int count_bonds2;
	int count_dis;

	// varies with time
	int time_nnn;
	double* time_dis1;
	double* time_dis1_max;
	double* time_dis1_min;
	double* time_count1;
	double* time_dis2;
	double* time_dis2_max;
	double* time_dis2_min;
	double* time_count2;
	double* time_accept;

	// mohan added 2017-03-19
	int np_tot; // total points
	double* pt_time; // time associated with pt
	double* coor_num;// coordination number associated with pt
};


#endif
