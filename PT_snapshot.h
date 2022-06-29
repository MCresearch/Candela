#ifndef PT_SS_H
#define PT_SS_H

#include "cellFile.h"
#include "HBs.h"
#include "input.h"

class PT_snapshot
{
	public:

	PT_snapshot();
	~PT_snapshot();

	void Routine();
	void print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index);
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

	int* snapshot_printed;
	int* snapshot_printed_pr;
	int* snapshot_printed_aft;
};

#endif