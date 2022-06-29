#include "cellFile.h"
#include "input.h"
#include "ili_3D.h"
#include "math.h"
#include "water.h"
#include "HBs.h"
#include "mdp.h"

ILI_3D::ILI_3D()
{
}

ILI_3D::~ILI_3D(){}

void ILI_3D::Routine()
{
	TITLE("ILI_3D","Routine");
	
	cout << "Compute the 3D Instantaneous Liquid Interface (ILI)" << endl;

	assert(INPUT.nx>0);
	assert(INPUT.ny>0);
	assert(INPUT.nz>0);

	// setup interface
	// zindex, interface, and gradient
	this->zindex = new int*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->zindex[ix] = new int[INPUT.ny];
	}
	this->interface = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->interface[ix] = new double[INPUT.ny];
	}
	this->gradient = new double**[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->gradient[ix] = new double*[INPUT.ny];
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			this->gradient[ix][iy] = new double[3]();
		}
	}


	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;


	// input ili file
	ifstream ifs(INPUT.ili_file.c_str());
	if(!ifs)
	{
		cout << "Cannot find the ILI file: " << INPUT.ili_file << endl;
		exit(0);
	}
	else
	{
		cout << "Open the ILI file: " << INPUT.ili_file << endl;
	}
	
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

//		cout << igeo << " use:" << cel.read_and_used << endl;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << " geo_target=" << INPUT.geo_target << endl;

		if(igeo==INPUT.geo_target)
		{
			plot(ifs, cel);
			break;
		}
	}	

	// delete arrays
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		delete[] interface[ix];
		delete[] zindex[ix];
	}
	delete[] interface;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			delete[] gradient[ix][iy];
		}
		delete[] gradient[ix];
	}
	delete[] gradient;

	return;
}

void ILI_3D::plot(ifstream &ifs, const Cell &cel)
{
	// type
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

	// water information
	Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;

    HBs::setup_water(cel, water);
	cout << "number of ions: " << Water::nions << endl;

	cout << "plot 3D ILI." << endl;

	cout << "Read the interface" << endl;
	MDP::read_ili(ifs, this->interface, this->gradient);

	ofstream ofs("ILI3D.xsf");
	ofs << "CRYSTAL" << endl;
	ofs << "PRIMVEC" << endl;
	ofs << INPUT.celldm1 << " 0 0" << endl;
	ofs << "0 " << INPUT.celldm2 << " 0" << endl;
	ofs << "0 0 " << INPUT.celldm3 << endl;
	ofs << "PRIMCOORD" << endl;
	ofs << INPUT.natom << " 1" << endl;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			if(INPUT.only_hydroxide==true and it==ito and water[ia].nH==1)
			{
				ofs << "S" << " " << cel.atom[it].pos[ia].x 
				<< " " << cel.atom[it].pos[ia].y
				<< " " << cel.atom[it].pos[ia].z << endl;	
			}
			else
			{
				ofs << cel.atom[it].id << " " << cel.atom[it].pos[ia].x 
				<< " " << cel.atom[it].pos[ia].y
				<< " " << cel.atom[it].pos[ia].z << endl;	
			}

			if(it==ito)
			{
				ofs_running << setw(10) << "water" << setw(10) << ia+1 
				<< setw(10) << water[ia].ndonate << setw(10) << water[ia].naccept << endl;
			}
		}
	} 
	ofs << "BEGIN_BLOCK_DATAGRID_3D" << endl;
	ofs << "3D_PWSCF" << endl;
	ofs << "DATAGRID_3D_UNKNOWN" << endl;
	ofs << INPUT.nx << " " << INPUT.ny << " " << INPUT.nz << endl;
	ofs << "0.000 0.000 0.000" << endl;
	ofs << INPUT.celldm1 << " 0 0" << endl;
	ofs << "0 " << INPUT.celldm2 << " 0" << endl;
	ofs << "0 0 " << INPUT.celldm3 << endl;
	

	assert(INPUT.nz>0);
	double dz = INPUT.celldm3/INPUT.nz;
	for(int iy=0; iy<INPUT.ny; ++iy)
	{
		for(int ix=0; ix<INPUT.nx; ++ix)
		{
			double zz = interface[ix][iy];
			int iz = (int)(zz/dz);
	//		cout << zz << " " << dz << " " << iz << endl;
			assert(iz<INPUT.nz and iz>=0);
			this->zindex[ix][iy]=iz;
		}
	}

	// z, y, x
	for(int iz=0; iz<INPUT.nz; ++iz)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			for(int ix=0; ix<INPUT.nx; ++ix)
			{
				if(zindex[ix][iy]==iz)
				{
					ofs << "1 ";
				}
				else ofs << "0 ";
			}
			ofs << endl;
		}
	}	

	ofs << "END_DATAGRID_3D" << endl;
	ofs << "END_BLOCK_DATAGRID_3D" << endl;
	

	ofs.close();
	delete[] water;
}
