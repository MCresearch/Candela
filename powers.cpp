#include "input.h"
#include "powers.h"

void PowerSpectra::Routine()
{
	TITLE("PowerSpectra","Routine");

	cal();

	return;
}

void PowerSpectra::cal() 
{
	TITLE("PowerSpectra","cal");

	// (1) open the file
	string name = INPUT.ps_in;
	
	ifstream ifs(name.c_str());
	if(!ifs)
	{
		cout << " Can't find file : " << name << endl;
		return;
	}
	cout << " Find file " << name << endl;

	// (2) read in VAF data
    int nv = INPUT.ps_nv;
	cout << " Number of points in this file : " << nv << endl;
	assert( nv > 0 );
	double* time = new double[nv];
	double* vaf = new double[nv];

	for(int i=0; i<nv; ++i)
	{
		ifs >> time[i] >> vaf[i];
		time[i] *= INPUT.ps_dt; // unit is fs
	}

	// (3) do FFT
	double dw = INPUT.ps_dw;
	double pi = 3.1415926535897;
	double two_pi = 2 * pi;
	int nw = INPUT.ps_nw;
	assert(dw>0.0); assert(nw>0);
	double* ps = new double[nw];

	double to_Herts = 2.99793e10; // equal 1 cm^-1

	double max_time = time[nv-1];
	assert(max_time > 0.0);
	cout << " maximal time (second) is " << max_time << endl;
	cout << " maximal frequency (cm^-1) could be " << two_pi/INPUT.ps_dt/to_Herts << endl;

	for(int i=0; i<nw; ++i) // i: points in powerspectra
	{
		// omega is in unit of cm^-1
		double omega = dw * i;
		// change the cm^-1 to Hertz
		omega *= to_Herts;
		ps[i] = 0.0;
		for(int j=0; j<nv; ++j) // j: points in VAF
		{
			// time is in unit of second
			// the frequency range is too narrow, unreasonable
			//ps[i] += vaf[j]*cos(two_pi * omega * time[j]);
			// should be the correct one
			ps[i] += vaf[j]*cos(omega * time[j]);
		}
	}
	
	// output the results
	ofstream ofs(INPUT.ps_out.c_str());
	if(ofs)
	{
		cout << " output results in " << INPUT.ps_out << endl;
	}
	for(int i=0; i<nw; ++i)
	{
		// i*dw is in unit of cm^-1
		ofs << setw(10) << i*dw << setw(15) << ps[i] << endl;
	}
	ofs.close();

	// clean up
	delete[] time;
	delete[] vaf;
	delete[] ps;
	return;
}


