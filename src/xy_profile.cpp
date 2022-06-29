#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "xy_profile.h"
#include "math.h"
#include "water.h"

XY_Profile::XY_Profile(){}
XY_Profile::~XY_Profile(){}

void XY_Profile::Routine()
{
	cout << "compute the XY profile"  << endl;
	compute();
	return;
}

void XY_Profile::compute()
{

	assert(INPUT.nx>0);
	assert(INPUT.ny>0);
	int nx = INPUT.nx;
	int ny = INPUT.ny;

	this->coord_xy = new double*[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->coord_xy[ix] = new double[ny]();
	}


	// print out the geometry
	this->count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;

		xy_coordinate(cel, igeo);
	}	

	ofstream ofs_result("xy_profile.dat");

	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			if(coord_xy[ix][iy]>0)
			{
				coord_xy[ix][iy]=-std::log(coord_xy[ix][iy]); 
			}
		}
	}

	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			ofs_result << coord_xy[ix][iy] << " "; 
		}
		ofs_result << endl;
	}
	ofs_result.close();



	// clean up
	for(int ix=0; ix<nx; ++ix)
	{
		delete[] coord_xy[ix];
	}
	delete[] coord_xy;


	return;
}


void XY_Profile::xy_coordinate(const Cell &cel, const int &igeo)
{
    int ito=-1;
    int ith=-1;
    int itc=-1;
    for(int it=0;it <INPUT.ntype; ++it)
    {
        if(cel.atom[it].id=="O") ito=it;
        else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
        else if(cel.atom[it].id=="C") itc=it;
    }
    if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
    if(INPUT.ntype==3){ assert(itc>=0); }

    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;

    HBs::setup_water(cel, water);

	if(INPUT.func_b==1)
	{
		for(int ia=0; ia<cel.atom[itc].na; ++ia)
		{
			double xxx = cel.atom[itc].pos[ia].x;
			double yyy = cel.atom[itc].pos[ia].y;
			double zzz = cel.atom[itc].pos[ia].z;
			while(xxx<0) xxx+=INPUT.celldm1;
			while(xxx>=INPUT.celldm1) xxx-=INPUT.celldm1;
			while(yyy<0) yyy+=INPUT.celldm2;
			while(yyy>=INPUT.celldm2) yyy-=INPUT.celldm2;
			while(zzz<0) zzz+=INPUT.celldm3;
			while(zzz>=INPUT.celldm3) zzz-=INPUT.celldm3;
			if( zzz > INPUT.z0 and zzz < INPUT.z1 )
			{
				int iX = INPUT.nx * xxx/INPUT.celldm1;
				int iY = INPUT.ny * yyy/INPUT.celldm2;
				assert(iX>=0 and iX<INPUT.nx);
				assert(iY>=0 and iY<INPUT.ny);
				coord_xy[iX][iY]+=1.0;
			}
		}
	}
	else if(INPUT.func_b==2)
	{
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if(water[ia].nH==1)
			{
				double xxx = cel.atom[ito].pos[ia].x;
				double yyy = cel.atom[ito].pos[ia].y;
				while(xxx<0) xxx+=INPUT.celldm1;
				while(xxx>=INPUT.celldm1) xxx-=INPUT.celldm1;
				while(yyy<0) yyy+=INPUT.celldm2;
				while(yyy>=INPUT.celldm2) yyy-=INPUT.celldm2;

				int iX = INPUT.nx * xxx/INPUT.celldm1;
				int iY = INPUT.ny * yyy/INPUT.celldm2;
				assert(iX>=0 and iX<INPUT.nx);
				assert(iY>=0 and iY<INPUT.ny);

				coord_xy[iX][iY]+=1.0;

				cout << "Hydroxide is on atom " << ia+1 << endl;
			}
		}
	}
	else if(INPUT.func_b==3)
	{
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if( cel.atom[ito].pos[ia].z > INPUT.z0 and cel.atom[ito].pos[ia].z < INPUT.z1 )
			{
				double xxx = cel.atom[ito].pos[ia].x;
				double yyy = cel.atom[ito].pos[ia].y;
				while(xxx<0) xxx+=INPUT.celldm1;
				while(xxx>=INPUT.celldm1) xxx-=INPUT.celldm1;
				while(yyy<0) yyy+=INPUT.celldm2;
				while(yyy>=INPUT.celldm2) yyy-=INPUT.celldm2;

				int iX = INPUT.nx * xxx/INPUT.celldm1;
				int iY = INPUT.ny * yyy/INPUT.celldm2;
				assert(iX>=0 and iX<INPUT.nx);
				assert(iY>=0 and iY<INPUT.ny);

				coord_xy[iX][iY]+=1.0;
			}
		}
	}
}
