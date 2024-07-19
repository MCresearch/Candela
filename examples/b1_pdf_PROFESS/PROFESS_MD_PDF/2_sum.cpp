#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cassert>
#include <ctime>

using namespace std;

int main()
{
	// sum up all the pair distribution function
	// from ion.start.dat.txt
	// to ion.end.dat.txt
	int start=0;
	int stop=0;
	int nmesh=0;
	int option=0;

	cout << " start, end, nmesh, option(1:pdf 2:ssf) : ";
    cin >> start >> stop >> nmesh >> option;
	
	assert(stop>=start);
	assert(nmesh>0);

	double* mesh = new double[nmesh];
	double* pdf = new double[nmesh];

	for(int i=0; i<nmesh; ++i)
	{
		mesh[i] = 0.0;
		pdf[i] = 0.0;
	}

	int count=0;
	for(int i=start; i<=stop; ++i)
	{
    	stringstream ss;
		if(option==1)
		ss << "ion." << i << ".dat.txt";
		else if(option==2)
		ss << "ion." << i << ".dat.ssf.txt";
		else
		{
			cout << "option should be 1 or 2." << endl;
			exit(0);
		}
		ifstream ifs(ss.str().c_str());

        if(!ifs)
		{
//        	cout << " Can't find the file: " << ss.str() << endl;
		}
	    else
		{
			if(option==2)
			{
				cout << " sum up ssf from file : " << ss.str() << endl; 
				for(int ir=0; ir<nmesh; ++ir)
				{
					double a,b;
					ifs >> a >> b;
					mesh[ir] = a;
					pdf[ir] += b;
				}
			}
			else if(option==1)
			{
				cout << " sum up pdf from file : " << ss.str() << endl; 
				for(int ir=0; ir<nmesh; ++ir)
				{
					double a,b,c;
					ifs >> a >> b ;
					mesh[ir] = a;
					pdf[ir] += b;
				}
			}
			++count;
		}

		ifs.close();
	}
	cout << "count=" << count << endl;

	string outfile;
	if(option==1) outfile="result-pdf.txt";
	else if(option==2) outfile="result-ssf.txt";

	ofstream ofs(outfile.c_str());
	assert(count!=0);
	for(int ir=0; ir<nmesh; ++ir)
	{
		ofs << mesh[ir] << " " << pdf[ir]/count << endl;
	}
	ofs.close();
	cout << " final result in : " << outfile << endl;

	delete[] mesh;
	delete[] pdf;

	return 0;
}
