#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "math.h"
#include "water.h"
#include "pdf.h" // to use compute_delta
#include "hyper.h"
#include "directional.h"

Directional::Directional()
{
}

Directional::~Directional(){}

void Directional::Routine()
{
	assert(INPUT.nx>0);
	assert(INPUT.ny>0);

	// initialize 2d array
	this->drx = 5.0/INPUT.nx; // unit is Angstroms
	this->dry = 40/INPUT.ny; // 180 degrees

	this->dist2D = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->dist2D[ix] = new double[INPUT.ny]();
	}

	this->count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		cout << "igeo=" << igeo << endl;

		compute_dist2D(cel);
		++count_geometry_number;
	}


	double sum=0.0;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			sum += dist2D[ix][iy];
		}
	}
	assert(sum>0.0);

	// print out the data
	ofstream ofs("directional.dat");
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			dist2D[ix][iy] = dist2D[ix][iy]/sum/drx/dry;
			ofs << dist2D[ix][iy] << " ";
		}
		ofs << endl;
	}
	ofs.close();


	// delete data
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		delete[] dist2D[ix];
	}
	delete[] dist2D;

	return;
}

void Directional::compute_dist2D(const Cell &cel)
{
	// get it index for each element;
	int ito=-1;
	int ith=-1;
	int itcl=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
	}
	if(INPUT.ntype>=2){ assert(ito>=0); assert(ith>=0);}

	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

	HBs::setup_water(cel, water);

	for(int i=0; i<cel.atom[ito].na; ++i)
	{
		for(int j=0; j<water[ito].ndonate; ++j)
		{
			double angle0 = water[ito].donate_angle[j];
			double dis = water[ito].donate_disO[j];
			locate(angle0,dis);
		}	
	}


	delete[] water;
}


void Directional::locate(const double &angle0, const double &dis)
{
	int indx = (double)dis/drx;
	int indy = (double)angle0/dry;

	if(indx<INPUT.nx and indy<INPUT.ny and indx>=0 and indy>=0)
	{
		dist2D[indx][indy] += 1;	
	}
	return;
}

