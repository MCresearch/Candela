#ifndef XSF_H
#define XSF_H

#include "cell.h"

class XSF 
{
	public: 
	
	XSF(){};
	~XSF(){};

	void Routine();
	
	private:

	void cal();

	void read_xsf(ifstream &ifs, const int &igeo);


};


#endif
