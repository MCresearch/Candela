#ifndef PT_SS2_H
#define PT_SS2_H

#include "cellFile.h"
#include "HBs.h"
#include "input.h"
#include "PT_snapshot.h"

class PT_snapshot2
{
public:
	PT_snapshot2();
	~PT_snapshot2();
	void print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index);
	void print_input_cp(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index);	
	void Routine();
	int* index;
	int* PT_start_snapshot_index;
	double* PT_start_snapshot_time;
	int* PT_end_snapshot_index;
	double* PT_end_snapshot_time;
	int* ion_index;
	int* nacc_start;
	int* ndon_start;
	int* nacc_end;
	int* ndon_end;
	int* H_index;
};

#endif