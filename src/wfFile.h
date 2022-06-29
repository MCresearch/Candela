#ifndef WFFILE_H
#define WFFILE_H

#include "rwstream.h"
#include "stdio.h"

#include "gfun.h"
#include "vec3.h"
#include "cell.h"

class WfFile  
{
	public: 
	
	WfFile();
	~WfFile();

	static bool ReadWf();

	private:

	static bool ReadWf_PWmat();
	private:

	static bool wf_open;
	static FILE *fp_kept;
};


#endif
/*
wf_in_type
wf_directory
