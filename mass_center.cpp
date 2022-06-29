#include "mass_center.h"
#include "cellFile.h"
#include "input.h"

mass_center::mass_center(){};

mass_center::~mass_center(){};

void mass_center::Routine()
{

	int count_geometry_number = 0;
	int na = INPUT.natom;
	pre_wpos = new Vector3<double>[na];
	wpos = new Vector3<double>[na];

	ofstream ofs("mass_center.txt");

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		double mx, my, mz;

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
		cout << "igeo=" << igeo << " time=" << cel.snapshot_time << endl;
		int ito=-1;
		int ith=-1;
		int itc=-1;

		//int na=0;
		for(int it=0;it <INPUT.ntype; ++it)
		{
			if(cel.atom[it].id=="O") ito=it;
			else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
			else if(cel.atom[it].id=="C") itc=it;
		//	na += cel.atom[it].na;
		}
		if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
		if(INPUT.ntype==3){ assert(itc>=0); }

		//int na = cel.atom[ito].na;
		
		//wpos0 = new Vector3<double>[na];
		assert(INPUT.celldm1>0);
		assert(INPUT.celldm2>0);
		assert(INPUT.celldm3>0);
		mx = 0;
		my = 0;
		mz = 0;
		double total_mass = 0;
		//double mass_moment = 0;
		if(count_geometry_number==1)
		{
			//this->start_time = cel.snapshot_time;
			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					wpos[iat] = cel.atom[it].pos[ia];
					pre_wpos[iat] = cel.atom[it].pos[ia];
					//wpos0[iat] = wpos[iat];
					++iat;
				}
			}
		}//end count_msd
		else
		{
			Vector3<double> move;
			int index0 = igeo-INPUT.geo_1;

			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					move.x = shortest(pre_wpos[iat].x, cel.atom[it].pos[ia].x, INPUT.celldm1);
					move.y = shortest(pre_wpos[iat].y, cel.atom[it].pos[ia].y, INPUT.celldm2);
					move.z = shortest(pre_wpos[iat].z, cel.atom[it].pos[ia].z, INPUT.celldm3);
					//wpos[iat] = wpos[iat] + move;  
					wpos[iat] = wpos[iat] - move;  // mohan update 2017-08-16
					pre_wpos[iat] = cel.atom[it].pos[ia];
					mx += wpos[iat].x*cel.atom[it].mass;
					my += wpos[iat].y*cel.atom[it].mass;
					mz += wpos[iat].z*cel.atom[it].mass;
					total_mass += cel.atom[it].mass;
					iat++;
				}
			}
		}
		ofs << cel.snapshot_index << " " << cel.snapshot_time << " " << mx/total_mass << " " << my/total_mass << " " << mz/total_mass << endl;
	}
	delete[] pre_wpos;
	delete[] wpos;
	ofs.close();
}