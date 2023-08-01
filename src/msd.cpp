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

	int na=0;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		na += cel.atom[it].na;
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
		pre_wpos = new Vector3<double>[na];
		wpos = new Vector3<double>[na];
		wpos0 = new Vector3<double>[na];
		allocate_water=true;
	}
	else if(INPUT.system=="oxygen" and allocate_water!=true)
	{	
		int na_o = cel.atom[ito].na;
		pre_wpos = new Vector3<double>[na_o];
		wpos = new Vector3<double>[na_o];
		wpos0 = new Vector3<double>[na_o];
		allocate_water=true;
	}
	else if(INPUT.system=="HPHT_water" and allocate_water!=true)
	{
		int na_o = cel.atom[ito].na;
		int na_h = cel.atom[ith].na;
		pre_wpos = new Vector3<double>[na_o];
		wpos = new Vector3<double>[na_o];
		wpos0 = new Vector3<double>[na_o];

		pre_wpos2 = new Vector3<double>[na_h];
		wpos2 = new Vector3<double>[na_h];
		wpos02 = new Vector3<double>[na_h];
		allocate_water=true;
	}
	if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		allocate_water=true;
	}
	else 
	{
		pre_wpos = new Vector3<double>[na];
		wpos = new Vector3<double>[na];
		wpos0 = new Vector3<double>[na];
	}

	Water *water;

	if (allocate_water)
	{
		water = new Water[cel.atom[ito].na];
		Water::nions=0;

		HBs::setup_water(cel, water);
	

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
		else if(INPUT.system=="oxygen")
		{
			if(count_msd==0)
			{
				this->start_time = cel.snapshot_time;
				for(int ia=0; ia<cel.atom[ito].na; ++ia)
				{
					if(water[ia].nH!=2) continue;
					//Vector3<double> mass_center = cel.atom[ito].pos[ia];

					double dx0 = shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[water[ia].indexH[0]].x, INPUT.celldm1);
					double dy0 = shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[water[ia].indexH[0]].y, INPUT.celldm2);
					double dz0 = shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[water[ia].indexH[0]].z, INPUT.celldm3);
					double dx1 = shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[water[ia].indexH[1]].x, INPUT.celldm1);
					double dy1 = shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[water[ia].indexH[1]].y, INPUT.celldm2);
					double dz1 = shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[water[ia].indexH[1]].z, INPUT.celldm3);
					double x0 = cel.atom[ito].pos[ia].x - dx0;	
					double y0 = cel.atom[ito].pos[ia].y - dy0;	
					double z0 = cel.atom[ito].pos[ia].z - dz0;	
					double x1 = cel.atom[ito].pos[ia].x - dx1;	
					double y1 = cel.atom[ito].pos[ia].y - dy1;	
					double z1 = cel.atom[ito].pos[ia].z - dz1;
					Vector3<double> H0 = Vector3<double>(x0,y0,z0);
					Vector3<double> H1 = Vector3<double>(x1,y1,z1);
					// mohan updated on 2018-12-24
					Vector3<double> mass_center = (cel.atom[ito].pos[ia]*cel.atom[ito].mass
						+H0*cel.atom[ith].mass+H0*cel.atom[ith].mass)
						/(cel.atom[ito].mass+2*cel.atom[ith].mass);

					wpos[ia] = mass_center;
					pre_wpos[ia] = wpos[ia];
					wpos0[ia] = wpos[ia];
				}
			}
			else
			{
				Vector3<double> move;
				for(int ia=0; ia<cel.atom[ito].na; ++ia)
				{
					if(water[ia].nH!=2) continue;
					//Vector3<double> mass_center = cel.atom[ito].pos[ia];

					double dx0 = shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[water[ia].indexH[0]].x, INPUT.celldm1);
					double dy0 = shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[water[ia].indexH[0]].y, INPUT.celldm2);
					double dz0 = shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[water[ia].indexH[0]].z, INPUT.celldm3);
					double dx1 = shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[water[ia].indexH[1]].x, INPUT.celldm1);
					double dy1 = shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[water[ia].indexH[1]].y, INPUT.celldm2);
					double dz1 = shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[water[ia].indexH[1]].z, INPUT.celldm3);
					double x0 = cel.atom[ito].pos[ia].x - dx0;	
					double y0 = cel.atom[ito].pos[ia].y - dy0;	
					double z0 = cel.atom[ito].pos[ia].z - dz0;	
					double x1 = cel.atom[ito].pos[ia].x - dx1;	
					double y1 = cel.atom[ito].pos[ia].y - dy1;	
					double z1 = cel.atom[ito].pos[ia].z - dz1;
					Vector3<double> H0 = Vector3<double>(x0,y0,z0);
					Vector3<double> H1 = Vector3<double>(x1,y1,z1);
					Vector3<double> mass_center = (cel.atom[ito].pos[ia]*cel.atom[ito].mass
						+H0*cel.atom[ith].mass+H0*cel.atom[ith].mass)
						/(cel.atom[ito].mass+2*cel.atom[ith].mass);

	//				Vector3<double> mass_center = 
	//					(cel.atom[ito].pos[ia]*16.0+
	//					cel.atom[ith].pos[water[ia].indexH[0]]*2.0+
	//					cel.atom[ith].pos[water[ia].indexH[1]]*2.0)/
	//					(16.0+2.0+2.0);
					move.x = shortest(pre_wpos[ia].x, mass_center.x, INPUT.celldm1);
					move.y = shortest(pre_wpos[ia].y, mass_center.y, INPUT.celldm2);
					move.z = shortest(pre_wpos[ia].z, mass_center.z, INPUT.celldm3);
					wpos[ia] = wpos[ia] + move; 
					pre_wpos[ia] = mass_center;
					double dx = wpos[ia].x-wpos0[ia].x;
					double dy = wpos[ia].y-wpos0[ia].y; 
					double dz = wpos[ia].z-wpos0[ia].z;
					this->msd[igeo-INPUT.geo_1] += dx*dx + dy*dy + dz*dz;  // unit is A^2
				}
				if(cel.snapshot_time>0.0)
				{
					ofs_msd << cel.snapshot_time << " " 
						<< msd[igeo-INPUT.geo_1]/cel.atom[ito].na << " " 
						<< msd[igeo-INPUT.geo_1]/6.0/(cel.snapshot_time-start_time)/cel.atom[ito].na << endl;
				}
			}
			++count_msd;
		}
	}
	else if (INPUT.system == "HPHT_water") // renxi added 20200912
	{
		if(count_msd==0)
		{
			this->start_time = cel.snapshot_time;

			for(int ia1=0; ia1<cel.atom[ito].na; ++ia1)
			{
				wpos[ia1] = cel.atom[ito].pos[ia1];
				pre_wpos[ia1] = cel.atom[ito].pos[ia1];
				wpos0[ia1] = wpos[ia1];
			}
			for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
			{
				wpos2[ia2] = cel.atom[ith].pos[ia2];
				pre_wpos2[ia2] = cel.atom[ith].pos[ia2];
				wpos02[ia2] = wpos2[ia2];
			}
		}//end count_msd
		else
		{
			Vector3<double> move;
		//	Vector3<double> move2;
			int index0 = igeo-INPUT.geo_1;
			for(int ia=0; ia<cel.atom[ito].na; ++ia)
			{
				move.x = shortest(pre_wpos[ia].x, cel.atom[ito].pos[ia].x, INPUT.celldm1);
				move.y = shortest(pre_wpos[ia].y, cel.atom[ito].pos[ia].y, INPUT.celldm2);
				move.z = shortest(pre_wpos[ia].z, cel.atom[ito].pos[ia].z, INPUT.celldm3);
				//wpos[iat] = wpos[iat] + move;  
				wpos[ia] = wpos[ia] - move;  // mohan update 2017-08-16
				pre_wpos[ia] = cel.atom[ito].pos[ia];
				double dx = wpos[ia].x-wpos0[ia].x;
				double dy = wpos[ia].y-wpos0[ia].y; 
				double dz = wpos[ia].z-wpos0[ia].z;
				double dxyz = dx*dx + dy*dy + dz*dz;
				//cout << iat << " " << dxyz << endl;
				this->msd[index0] += dxyz;  // unit is A^2
			}
			for(int ia=0; ia<cel.atom[ith].na; ++ia)
			{
				move.x = shortest(pre_wpos2[ia].x, cel.atom[ith].pos[ia].x, INPUT.celldm1);
				move.y = shortest(pre_wpos2[ia].y, cel.atom[ith].pos[ia].y, INPUT.celldm2);
				move.z = shortest(pre_wpos2[ia].z, cel.atom[ith].pos[ia].z, INPUT.celldm3);
				//wpos[iat] = wpos[iat] + move;  
				wpos2[ia] = wpos2[ia] - move;  // mohan update 2017-08-16
				pre_wpos2[ia] = cel.atom[ith].pos[ia];
				double dx = wpos2[ia].x-wpos02[ia].x;
				double dy = wpos2[ia].y-wpos02[ia].y; 
				double dz = wpos2[ia].z-wpos02[ia].z;
				double dxyz = dx*dx + dy*dy + dz*dz;
				//cout << iat << " " << dxyz << endl;
				this->msd2[index0] += dxyz;  // unit is A^2
			}
		}//end if
		if(cel.snapshot_time>0.0)
		{
			ofs_msd << cel.snapshot_time << " " 
					<< msd[igeo-INPUT.geo_1]/cel.atom[ito].na << " " 
					<< msd2[igeo-INPUT.geo_1]/cel.atom[ith].na << endl;
		}
		count_msd++;
	}
	else
	{
		if(count_msd==0)
		{
			this->start_time = cel.snapshot_time;
			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					wpos[iat] = cel.atom[it].pos[ia];
					pre_wpos[iat] = cel.atom[it].pos[ia];
					wpos0[iat] = wpos[iat];
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
					double dx = wpos[iat].x-wpos0[iat].x;
					double dy = wpos[iat].y-wpos0[iat].y; 
					double dz = wpos[iat].z-wpos0[iat].z;
					double dxyz = dx*dx + dy*dy + dz*dz;
					//cout << iat << " " << dxyz << endl;
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
