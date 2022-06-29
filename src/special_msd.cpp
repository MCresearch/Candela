#include "cellFile.h"
#include "input.h"
#include "special_msd.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"

special_MSD::special_MSD(){}
special_MSD::~special_MSD(){}

void special_MSD::Routine()
{
	TITLE("special_MSD","Routine");
	
	cout << "Compute the Mean Square Displacements of a specific atom ia_select" << endl;
	assert(INPUT.ia_select>=0);
	this->nmsd = INPUT.geo_2-INPUT.geo_1+1;	
	assert(nmsd>0);
	this->msd = new double[nmsd]();
	this->msd2 = new double[nmsd]();	
	this->count_msd = 0;
	this->start_time = 0.0;

	allocate_water=false;

	// setup the proton transfer data
	PT.setup_PT();

	ofstream ofs_msd("MSD.dat");

	//**************************
	// BEGIN CALCULATING DATA
	//**************************
	this->count_geometry_number=0;
	int ipt=0;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

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

		int na = cel.atom[ito].na;

		assert(INPUT.celldm1>0);
		assert(INPUT.celldm2>0);
		assert(INPUT.celldm3>0);
		if(INPUT.system=="water" and allocate_water==false)
		{
			pre_wpos = new Vector3<double>[na];
			wpos = new Vector3<double>[na];
			wpos0 = new Vector3<double>[na];
			allocate_water=true;
		}
		//else
		//{
		//	cout << "Something wrong, please check." << endl;
		//	exit(0);
		//}
//		Water *water = new Water[cel.atom[ito].na];
//		Water::nions=0;
//		HBs::setup_water(cel, water);
		int ia_select = INPUT.ia_select;
		if(count_msd==0)
		{
			this->start_time = cel.snapshot_time;
			for(int ia=0; ia<na; ia++)
			{
				wpos[ia] = cel.atom[ito].pos[ia];
				pre_wpos[ia] = cel.atom[ito].pos[ia];
				wpos0[ia] = wpos[ia];
			//	cout << wpos[0].x << " " << wpos[0].y << " " << wpos[0].z << endl;
			}
			
		}//end count_msd
		else
		{
			Vector3<double> move;
			int index0 = igeo-INPUT.geo_1;
			ofs_msd << cel.snapshot_time << " ";
			for(int ia=0; ia<na; ia++)
			{
				move.x = shortest(pre_wpos[ia].x, cel.atom[ito].pos[ia].x, INPUT.celldm1);
				move.y = shortest(pre_wpos[ia].y, cel.atom[ito].pos[ia].y, INPUT.celldm2);
				move.z = shortest(pre_wpos[ia].z, cel.atom[ito].pos[ia].z, INPUT.celldm3);
				//cout << cel.atom[ito].pos[ia_select].x << " " << cel.atom[ito].pos[ia_select].y << " " << cel.atom[ito].pos[ia_select].z << endl;
				//cout << pre_wpos[0].x << " " << pre_wpos[0].y << " " << pre_wpos[0].z << endl;
				//cout << move.x << " " << move.y << " " << move.z << endl;
				//wpos[iat] = wpos[iat] + move;  
				wpos[ia] = wpos[ia] - move;  // mohan update 2017-08-16
				pre_wpos[ia] = cel.atom[ito].pos[ia];
				double dx = wpos[ia].x-wpos0[ia].x;
				double dy = wpos[ia].y-wpos0[ia].y; 
				double dz = wpos[ia].z-wpos0[ia].z;
				double dxyz = dx*dx + dy*dy + dz*dz;
				//cout << iat << " " << dxyz << endl;
				this->msd[index0] = dxyz;  // unit is A^2
				ofs_msd << dxyz << " ";
			}
			ofs_msd << endl;
		}
		//if(cel.snapshot_time>0.0)
		//	{
				 
					
		//	}
		count_msd++;
	} //end geo
	ofs_msd.close();
}
