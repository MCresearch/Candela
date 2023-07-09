#include "cellFile.h"
#include "input.h"
#include "gfun.h"
#include "string.h"

bool CellFile::ReadGeometry_ABACUS_old( Cell &cel, ifstream &ifs )
{
	TITLE("CellFile","ReadGeometry_ABACUS_old");
	const int ntype = INPUT.ntype;
	bool restart = true;

	string useless;
	double lat_const;
	ifs >> useless >> useless;
	ifs >> useless >> lat_const;
	ifs >> useless;
	ifs >> cel.a1.x >> cel.a1.y >> cel.a1.z;
	ifs >> cel.a2.x >> cel.a2.y >> cel.a2.z;
	ifs >> cel.a3.x >> cel.a3.y >> cel.a3.z;

	cel.a1 *= lat_const*BOHR;
	cel.a2 *= lat_const*BOHR;
	cel.a3 *= lat_const*BOHR;

	INPUT.celldm1 = cel.a1.norm();
	INPUT.celldm2 = cel.a2.norm();
	INPUT.celldm3 = cel.a3.norm();
	cel.cal_volume();


	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt; 

	char keyword[20];
	ifs >> keyword;
	if (strcmp(keyword, "VIRIAL") == 0) 
	{
		for (int i=0; i<15; i++) ifs >> useless;
	}
	else if (strcmp(keyword, "INDEX") == 0)
	{
		for (int i=0; i<4; i++) ifs >> useless;
	}

	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		cel.nat+=cel.atom[it].na;
	}

	for(int it=0; it<ntype; ++it)
	{
		cel.atom[it].read_pos_5(ifs, cel.a1, cel.a2, cel.a3, lat_const);

		for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
		{
			cel.cartesian2direct(it, ia2);
		}
	}
	cel.atom_mass();
	return true;
}

