#include "cellFile.h"
#include "input.h"
#include "gfun.h"
//qianrui
bool CellFile::ReadGeometry_PWmat( Cell &cel, ifstream &ifs )
{
	TITLE("CellFile", "ReadGeometry_PWmat");
	const int ntype = INPUT.ntype;
	cel.init_cel(INPUT);
	const double celldm1 =INPUT.celldm1;
	const double celldm2 =INPUT.celldm2;
	const double celldm3 =INPUT.celldm3;
	string useless;
	string txt;
	while(ifs >> useless)
	{
		getline(ifs,txt);
		//cout << txt << endl;
		if(useless=="Lattice")
		{
			break;
		}

	}
	//read lattice
	ifs >> cel.a1.x >> cel.a1.y >> cel.a1.z;
	getline(ifs,txt);
	ifs >> cel.a2.x >> cel.a2.y >> cel.a2.z;
	getline(ifs,txt);
	ifs >> cel.a3.x >> cel.a3.y >> cel.a3.z;
	getline(ifs,txt);
	
	INPUT.celldm1 = sqrt(cel.a1.x*cel.a1.x + cel.a1.y*cel.a1.y + cel.a1.z*cel.a1.z);
    INPUT.celldm2 = sqrt(cel.a2.x*cel.a2.x + cel.a2.y*cel.a2.y + cel.a2.z*cel.a2.z);
    INPUT.celldm3 = sqrt(cel.a3.x*cel.a3.x + cel.a3.y*cel.a3.y + cel.a3.z*cel.a3.z);

	assert(INPUT.celldm1>0);
	assert(INPUT.celldm2>0);
	assert(INPUT.celldm3>0);
	//calculate the volume of the cell
	cel.cal_volume();
	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt; 

	delete[] cel.atom;
	cel.atom = new Atoms[ntype];

	//locate atomic positions
	while(ifs>>useless)
	{
		getline(ifs,txt);
		if(useless=="Position")
		{
			break;
		}
	}
	if(ifs.eof())
	{
		cout<<"Geo isn't enough!"<<endl;
		return false;
	}

	cel.atom_mass();



	double tmpx,tmpy,tmpz;
	for(int it=0; it<ntype; ++it)
	{
		cel.atom[it].pos = new Vector3<double>[cel.atom[it].na];
		for(int ia=0;ia<cel.atom[it].na;ia++)
		{
			ifs >> useless >> tmpx >> tmpy >> tmpz >> useless >> useless >> useless;
			cel.atom[it].pos[ia].x=tmpx*celldm1;
			cel.atom[it].pos[ia].y=tmpy*celldm2;
			cel.atom[it].pos[ia].z=tmpz*celldm3;
		}
	}
	


	return true;
}
