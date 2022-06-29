#ifndef RWSTREAM_H
#define RWSTREAM_H

#include <stdio.h>
#include <complex>
class RWSTREAM
{
	public:
		RWSTREAM(){};
		RWSTREAM(FILE* &ptr);
		~RWSTREAM(){};
		void setptr(FILE* &ptr);
		FILE* fileptr;
};
//For template operator, you had better write defination in .h file!!
template<class T>
RWSTREAM& operator<<(RWSTREAM& wstream,T& data)
{
    int size=sizeof(data);
    fwrite(&data,size,1,wstream.fileptr);
    return wstream;
}
template<class T>
RWSTREAM& operator>>(RWSTREAM& rstream,T& data)
{
	int size=sizeof(data);
    fread(&data,size,1,rstream.fileptr);
    return rstream;
}

#endif
