#include "cellFile.h"
#include "input.h"
#include "xsf.h"
#include "math.h"
#include "HBs.h"

void XSF::Routine()
{
	TITLE("XSF","Routine");

	cal();

	return;
}


void XSF::cal()
{
	TITLE("XSF","cal");

	cout << "Compute the XSF related properties" << endl;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		cout << "igeo=" << igeo << endl;
		stringstream ss;
		ss << "KS_" << igeo <<  ".xsf";
		
		ifstream ifs(ss.str().c_str());
		if(!ifs)
		{
			cout << "could not find file KS_" << igeo << ".xsf" << endl;
			continue;
		}

		read_xsf(ifs, igeo);
		

		ifs.close();
	}

	return;
}

void XSF::read_xsf(ifstream &ifs, const int &igeo)
{
	string name;
	ifs >> name;
//	cout << name << endl;

	ifs >> name;
//	cout << name << endl;

	// cell in Angstroms
	double e11, e12, e13;
	double e21, e22, e23;
	double e31, e32, e33;

	ifs >> e11 >> e12 >> e13;
	ifs >> e21 >> e22 >> e23;
	ifs >> e31 >> e32 >> e33;

//	cout << "cell lengths: " << e11 << " " << e22 << " " << e33 << endl;  

	double volume = e11 * e22 * e33;
//	cout << "volume is " << volume << endl;

	ifs >> name;
//	cout << name << endl;
	
	int natom;
	ifs >> natom;
	cout << "number of atoms: " << natom << endl;

	int useless;
	ifs >> useless;

	string* id = new string[natom];
	double* cx = new double[natom];
	double* cy = new double[natom];
	double* cz = new double[natom];

	for(int i=0; i<natom; ++i)
	{
		ifs >> id[i] >> cx[i] >> cy[i] >> cz[i];
//		cout << id[i] << endl;
	}

	ifs >> name;
//	cout << name << endl;

	ifs >> name;
//	cout << name << endl;

	ifs >> name;
//	cout << name << endl;

	int nx, ny, nz;
	ifs >> nx >> ny >> nz;
	cout << "dims: " << nx << " " << ny << " " << nz << endl;
	int nxyz = nx*ny*nz;

	double number;
	for(int i=0; i<12; ++i) ifs >> number;

	// initialize the wave functions
	double*** wf = new double**[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		wf[ix] = new double*[ny];
		for(int iy=0; iy<ny; ++iy)
		{
			wf[ix][iy] = new double[nz];
		}
	}

	// make sure the OO atom is within the cell
	while(cx[0]<0.0) cx[0]+=e11;
	while(cy[0]<0.0) cy[0]+=e22;
	while(cz[0]<0.0) cz[0]+=e33;
	while(cx[0]>=e11) cx[0]-=e11;
	while(cy[0]>=e22) cx[0]-=e22;
	while(cz[0]>=e33) cx[0]-=e33;

	// open file
	stringstream ss;
	ss << "rho1d_" << igeo << ".dat";
	ofstream ofs(ss.str().c_str());

	double dr = INPUT.dr;
	assert(dr>0.0);
	double rcut = INPUT.rcut;
	assert(rcut>0.0);

	// number of radial mesh grids.
	int nmesh = int(rcut / dr) +  1;
	
	// radial distribution function.
	double* rho1d = new double[nmesh];

	ofs_running << " dr = " << dr << " Angstrom" << endl;
	ofs_running << " rcut = " << rcut << " Angstrom" << endl;
	ofs_running << " nmesh = " << nmesh << endl;

	double sum=0.0;
	double sum1=0.0;
	double sum2=0.0;
	double sum3=0.0;
	double r1 = 2.75;
	double r2 = 4.40;

	for(int iz=0; iz<nz; ++iz)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			for(int ix=0; ix<nx; ++ix)
			{
				ifs >> wf[ix][iy][iz];
				double w2 = pow(wf[ix][iy][iz],2);
	
				double xx = e11*(double)ix/nx;
				double yy = e22*(double)iy/ny;
				double zz = e33*(double)iz/nz;

				double dx = shortest(cx[0], xx, e11);	
				double dy = shortest(cy[0], yy, e22);	
				double dz = shortest(cz[0], zz, e33);	

				double dis = sqrt(dx*dx+dy*dy+dz*dz);
			
				if(dis < rcut)
				{
					int index = dis/dr;
					rho1d[index] += w2; 
				}

				if(dis < r1)
				{
					sum1 += w2;
				}
				else if(dis >= r1 and dis < r2)
				{
					sum2 += w2;
				}
				else if(dis >= r2)
				{
					sum3 += w2;
				}

				sum += w2;
			}
		}
	}

	for(int ir=0; ir<nmesh; ++ir)
	{
		ofs << (ir+0.5)*dr << " " << rho1d[ir]/(double)nxyz << endl; 
	}
	ofs.close();
	
	sum /= (double)nxyz;
	cout << "<psi|psi> = " << sum << endl;

	sum1 /= (double)nxyz;
	sum2 /= (double)nxyz;
	sum3 /= (double)nxyz;
	ofs_running << "decomposed_sum " << sum1 << " " << sum2 << " " << sum3 << " " << endl;

	// clean up
	delete[] id;
	delete[] cx;
	delete[] cy;
	delete[] cz;
	delete[] rho1d;

	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			delete[] wf[ix][iy];
		}
		delete[] wf[ix];
	}
	delete[] wf;

	return;
}
