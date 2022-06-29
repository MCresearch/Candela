#ifndef HBS_NEAR_PT_H
#define HBS_NEAR_PT_H

class HBs_near_PT
{
public:
	HBs_near_PT();
	~HBs_near_PT();

	void Routine();
	static int round(const double &r);
	double* hbs_bf_pt;
	double* hbs_aft_pt;
	int* count_bf_pt;
	int* count_aft_pt;
	int* snapshot_pt;
	double* snapshot_time_pt;
	int* iindex;
	int* iindex_p; 
	bool* return_jump;
};

#endif