#ifndef binfstream_H
#define binfstream_H
#include "malloc.h"
#include <stdio.h>
#include <cstdlib>
#include <complex>
#include <iostream>
using namespace std;
//A class to read or write binary files.
//By qianrui 2020-1-6
class binfstream
{
	public:
		binfstream(){
			fileptr=nullptr;
		};
		binfstream(const string,const char*);
		binfstream(FILE* &ptr);
		~binfstream(){};
		void setptr(FILE* &ptr);
		FILE* fileptr;
		void close();
		void open(const string,const char*);
		bool operator!() const;
		operator bool() const;
};
//For template operator, you had better write defination in .h file!!
//for variation or array
template<class T>
binfstream& operator<<(binfstream& wstream,const T data)
{
    int size=sizeof(T);
	int n=sizeof(data)/sizeof(T);
    fwrite(&data,size,n,wstream.fileptr);
    return wstream;
}
/*//for dynamic memory
//malloc_usable_size has problem!
template<class T>
binfstream& operator<<(binfstream& wstream,T* &data)
{
    int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
    fwrite(data,size,n,wstream.fileptr);
    return wstream;
}*/
template<class T>
binfstream& operator>>(binfstream& rstream,T& data)
{
	int size=sizeof(T);
	int n=sizeof(data)/sizeof(T);
	size_t ch;
    ch=fread(&data,size,n,rstream.fileptr);
#ifdef __JUG
	if(ch<n)
	{
		cout<<"Error in binfstream: Some data didn't be read."<<endl;
    	exit(0);
	}
#endif
	return rstream;
}
/*//for dynamic memory
template<class T>
binfstream& operator>>(binfstream& rstream,T* &data)
{
	int size=sizeof(T);
	int n=malloc_usable_size(data)/sizeof(T);
	cout<<malloc_usable_size(data)<<' '<<sizeof(T)<<' '<<n<<endl;
	size_t ch;
    ch=fread(data,size,n,rstream.fileptr);
	if(ch<n) 
	{
		cout<<"Error in binfstream: Some dynamic memory didn't be read."<<endl;
		exit(0);
    }
	return rstream;
}*/
template<class T>
void rwread(binfstream& rstream,T* &data,int n)
{
	int size=sizeof(T);
	size_t ch;
    ch=fread(data,size,n,rstream.fileptr);
#ifdef __JUG
    if(ch<n) 
    {
        cout<<"Error in binfstream: Some dynamic memory didn't be read."<<endl;
        exit(0);
    }
#endif
    return;
}	
#endif

