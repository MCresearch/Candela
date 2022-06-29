#ifndef WATERWIRE2_H
#define WATERWIRE2_H

#include "cell.h"
#include "water.h"
#include "mj.h"

class Wire
{
	public:
	
	Wire();
	~Wire();

	int o1;
	int o2;
	int o3;
	int h1;
	int h2;
	static int n;
	bool f_active;
};



// Properties of Waterwire 
class Waterwire2
{
	public: 
	
	Waterwire2();
	~Waterwire2();

	void Routine();

	private:

	int count_geometry_number;
	ProtonTransfer PT;

	void search_compression(const Cell &cel, const int &igeo);
	void double_jump(const Cell &cel, Water *water, const int &ito, const int &ith);
	void read_wire();
	void all_ions_involved_water(const Cell &cel, Water *water, const int &ito, const int &ith);
	void all_water(const Cell &cel, Water *water, const int &ito, const int &ith);

	double** coord_xy;
	double** freeE;
	ifstream ifs_wire;

	Wire* ww;
	int ss_index;
	double ss_time;

	double* dis_o123;
	double avg123;
	int count123;
	
	double avg_d;
	int count_d;
};


#endif
