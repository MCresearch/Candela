#include "cellFile.h"
#include "input.h"
#include "msd.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"

MSD::MSD(){}
MSD::~MSD(){}

void MSD::Routine()
{
	TITLE("MSD","Routine");
	
	cout << "Compute the Mean Square Displacements" << endl;

	this->nmsd = INPUT.geo_2-INPUT.geo_1+1;	
	assert(nmsd>0);
	assert(INPUT.msd_dt>0.0);

	this->msd = new double[nmsd]();
	this->msd2 = new double[nmsd]();	
	this->count_msd = 0;
	this->start_time = 0.0;

	for (int imsd = 0; imsd < nmsd; imsd++)
	{
		this->msd[imsd] = 0;
		this->msd2[imsd] = 0;
	}

		allocate_water = false;

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

		compute_msd(cel, igeo, ofs_msd);

		cel.clean();
	}	

	if(allocate_water)
	{
		delete[] pre_wpos;
		delete[] wpos;
		delete[] wpos0;
	}

	delete[] msd;
	ofs_msd.close();

	return;
}


void MSD::compute_msd(const Cell &cel, const int &igeo, ofstream &ofs_msd)
{
	// get ito, ith, and itc.
	int ito=-1;
	int ith=-1;
	int itc=-1;
	int it_select = -1;

	int na=0;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		na += cel.atom[it].na;

		if (INPUT.ele_select != "none")
		{
			if (cel.atom[it].id == INPUT.ele_select)
			{
				it_select = it;
			}
		}
	}
//	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
//	if(INPUT.ntype==3){ assert(itc>=0); }

	assert(INPUT.celldm1>0);
	assert(INPUT.celldm2>0);
	assert(INPUT.celldm3>0);

	int its = ito;
	if(INPUT.func==2)
	{
		its = ith;
	}


	if(INPUT.system=="water" and allocate_water!=true)
	{
		if (count_msd == 0)
		{
			pre_wpos = new Vector3<double>[na];
			wpos = new Vector3<double>[na];
			wpos0 = new Vector3<double>[na];
		}
		allocate_water=true;
	}
	if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		allocate_water=true;
	}
	else 
	{
		if (count_msd == 0 and it_select == -1)
		{
			pre_wpos = new Vector3<double>[na];
			wpos = new Vector3<double>[na];
			wpos0 = new Vector3<double>[na];
		}
		else if (count_msd == 0 and it_select != -1)
		{
			assert(it_select >= 0);
			pre_wpos = new Vector3<double>[cel.atom[it_select].na];
			wpos = new Vector3<double>[cel.atom[it_select].na];
			wpos0 = new Vector3<double>[cel.atom[it_select].na];
		}
	}

	Water *water;

	if (allocate_water)
	{
		water = new Water[cel.atom[ito].na];
		Water::nions=0;

		HBs::setup_water(cel, water);
	}

	if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if( INPUT.system=="hydronium" and water[ia].nH!=3) continue;
			if( INPUT.system=="hydroxide" and water[ia].nH!=1) continue;

			if(Water::nions==1)
			{
				if(count_msd==0)
				{
					ion_pos.x = cel.atom[ito].pos[ia].x;
					ion_pos.y = cel.atom[ito].pos[ia].y;
					ion_pos.z = cel.atom[ito].pos[ia].z;
					pre_pos = ion_pos;
					ion_pos0 = ion_pos;
					cout << setw(10) << cel.snapshot_time 
						<< setw(10) << ion_pos.x 
						<< setw(10) << ion_pos.y 
						<< setw(10) << ion_pos.z << endl;
				}
				else
				{
					Vector3<double> move;
					move.x = shortest(pre_pos.x, cel.atom[ito].pos[ia].x, INPUT.celldm1);
					move.y = shortest(pre_pos.y, cel.atom[ito].pos[ia].y, INPUT.celldm2);
					move.z = shortest(pre_pos.z, cel.atom[ito].pos[ia].z, INPUT.celldm3);
					ion_pos = ion_pos + move; 
					pre_pos = cel.atom[ito].pos[ia];
					double dx = ion_pos.x-ion_pos0.x;
					double dy = ion_pos.y-ion_pos0.y; // fix a bug 2017-02-08 by mohan
					double dz = ion_pos.z-ion_pos0.z;
					this->msd[igeo-INPUT.geo_1] = dx*dx + dy*dy + dz*dz;  // unit is A^2
					ofs_msd << setprecision(12) << cel.snapshot_time << " " << msd[igeo-INPUT.geo_1] << endl;

				}
				++count_msd;
			}
		}
	}
	else
	{
		if(count_msd==0)
		{
			this->start_time = cel.snapshot_time;
			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				if (it_select != -1 and it != it_select)
				{
					continue;
				}
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					this->wpos[iat] = cel.atom[it].pos[ia];
					this->pre_wpos[iat] = cel.atom[it].pos[ia];
					this->wpos0[iat] = wpos[iat];
					// cout << this->wpos0[iat].x << " " << wpos0[iat].y << " " << wpos0[iat].z << endl;
					++iat;
				}
			}
		} // end count_msd
		else
		{
			Vector3<double> move;
			int index0 = igeo-INPUT.geo_1;

			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				if (it_select != -1 and it != it_select)
				{
					continue;
				}
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					move.x = shortest(pre_wpos[iat].x, cel.atom[it].pos[ia].x, INPUT.celldm1);
					move.y = shortest(pre_wpos[iat].y, cel.atom[it].pos[ia].y, INPUT.celldm2);
					move.z = shortest(pre_wpos[iat].z, cel.atom[it].pos[ia].z, INPUT.celldm3);
					// cout << this->pre_wpos[iat].x << " " << this->pre_wpos[iat].y << " " << this->pre_wpos[iat].z << endl;
					// cout << cel.atom[it].pos[ia].x << " " << cel.atom[it].pos[ia].y << " " << cel.atom[it].pos[ia].z << endl;
					// cout << this->wpos0[iat].x << " " << this->wpos0[iat].y << " " << this->wpos0[iat].z << endl;
					// wpos[iat] = wpos[iat] + move;
					wpos[iat] = wpos[iat] - move;  // mohan update 2017-08-16
					pre_wpos[iat] = cel.atom[it].pos[ia];
					double dx = wpos[iat].x-wpos0[iat].x;
					double dy = wpos[iat].y-wpos0[iat].y; 
					double dz = wpos[iat].z-wpos0[iat].z;
					double dxyz = dx*dx + dy*dy + dz*dz;
					//cout << iat << " " << move.x << " " << move.y << " " << move.y << " " << dxyz << endl;
					this->msd[index0] += dxyz;  // unit is A^2
					++iat;
				}
			}
			if(cel.snapshot_time>0.0)
			{
				ofs_msd << cel.snapshot_time << " " 
					<< msd[index0]/iat << " " 
					<< msd[index0]/6.0/(cel.snapshot_time-start_time)/iat << endl;
			}
		}
		++count_msd;
	}//end system


	if (allocate_water)
	{
		delete[] water;
	}
	return;
}
