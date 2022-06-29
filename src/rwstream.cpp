#include <stdio.h>
#include <complex.h>
#include "rwstream.h"

RWSTREAM::RWSTREAM(FILE* &ptr)
{
	fileptr=ptr;
}

void RWSTREAM:: setptr(FILE* &ptr)
{
	fileptr=ptr;
	return;
}

