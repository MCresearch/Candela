#include <stdio.h>
#include "binfstream.h"

binfstream::binfstream(FILE* &ptr)
{
	fileptr=ptr;
}

binfstream::binfstream(const string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

void binfstream:: setptr(FILE* &ptr)
{
	fileptr=ptr;
	return;
}

void binfstream::open(const string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

void binfstream:: close()
{
	fclose(fileptr);
	return;
}

bool binfstream::operator!() const
{
	if (fileptr==nullptr)
		return true;
	else
		return false;
}
binfstream::operator bool() const
{
	if (fileptr==nullptr)
		return false;
	else
		return true;
}
