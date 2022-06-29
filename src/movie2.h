#ifndef MOVIE2_H
#define MOVIE2_H
#include "HBs.h"
#include "cellFile.h"
class movie2
{
public:
	~movie2(){};
	movie2(){};

	void Routine();
	void snapshot(ofstream &ofs,Cell &cel,int &igeo, int* &Ochain, int* &Hindex, bool &init, int &nO, int &nH);
};

#endif
