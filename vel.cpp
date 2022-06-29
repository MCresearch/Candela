#include "cellFile.h"
#include "input.h"
#include "vel.h"

void Vel::Routine()
{
	TITLE("Vel","Routine");
	// cel_in : input geometry file
	CellFile cel_in;

	cal( cel_in );

	return;
}


void Vel::cal( Cell &cel )
{
	TITLE("Vel","cal");

	string line1;
	string line2;
	string line3;
	string line4;

	double* vx = new double[INPUT.natom];
	double* vy = new double[INPUT.natom];
	double* vz = new double[INPUT.natom];
	double* vv = new double[INPUT.natom];

	for(int i=0; i<INPUT.natom; ++i)
	{
		vx[i] = 0.0;
		vy[i] = 0.0;
		vz[i] = 0.0;
		vv[i] = 0.0;
	}

	// read in the velocities.
	ifstream ifs(INPUT.vel_in.c_str());
	if(!ifs)
	{
		cout << " Can't find the file: " << INPUT.vel_in << endl;
		exit(0);
	}

	getline(ifs,line1);
	getline(ifs,line2);
	getline(ifs,line3);
	getline(ifs,line4);

	for(int i=0; i<INPUT.natom; ++i)
	{
		ifs >> vx[i] >> vy[i] >> vz[i];
		vv[i] = sqrt( vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );
	}


	// reorder the velocity
	double sx, sy, sz, sv;
	for(int i=0; i<INPUT.natom-1; ++i)
	{
		for(int j=i+1; j<INPUT.natom; ++j)
		{
			if(vv[i] > vv[j])
			{
				sx = vx[i]; sy = vy[i]; sz = vz[i]; sv = vv[i]; 
                vx[i] = vx[j]; vy[i] = vy[j]; vz[i] = vz[j];  vv[i] = vv[j];
				vx[j] = sx; vy[j] = sy; vz[j] = sz; vv[j] = sv;
			}
		}
	}


	cout << " The max velocity is " << vv[INPUT.natom-1] << endl;
	cout << " The min velocity is " << vv[0] << endl;

	assert( INPUT.ndv > 0 );
	int* count = new int[INPUT.ndv];
	double dv = (vv[INPUT.natom-1]-vv[0])/INPUT.ndv;
	for(int id=0; id<INPUT.ndv; ++id)
	{
		count[id]=0;
		double bound = dv*id;
		for(int i=0; i<INPUT.natom; ++i)
		{
			if(vv[i]>=bound && vv[i]<bound+dv)
			{
				count[id]++;
			}
		}
	}


    // output the final pair distribution function
	ofstream ofs(INPUT.vel_out.c_str());

//	// the search for last element is not complete,
//	// so we abandom the last point.
//	for(int i=0; i<INPUT.natom; ++i)
//	{
//		ofs << vv[i] << " " << vx[i] << " " << vy[i] << " " << vz[i] << endl;
//	}
	for(int id=0; id<INPUT.ndv; ++id)
	{
		ofs << id+1 << " " << count[id] << endl;
	}

	ifs.close();
	ofs.close();

	delete[] vx;
	delete[] vy;
	delete[] vz;
	delete[] count;
	return;
}


