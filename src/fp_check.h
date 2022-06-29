#ifndef FP_CHECK_H
#define FP_CHECK_H
#include "cellFile.h"
#include "HBs.h"
#include "input.h"

class fp_check
{
public:
	fp_check();
	~fp_check();

	void Routine();
	void print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index);
	void print_input_cp(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index);
};

#endif