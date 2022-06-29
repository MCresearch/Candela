#include "movie.h"
#include "input.h"
#include "atoms.h"
#include "HBs.h"

Movie::Movie() 
{
}

Movie::~Movie() 
{
}

void Movie::Routine()
{
	TITLE("Movie","Routine");
	
	ofs_running << "Produce a movie for MD trajectory" << endl;

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	// output data
	ofstream ofs(INPUT.geo_out.c_str());

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//ofs_running << "geometry " << igeo%INPUT.geo_interval << endl;
		
		if(igeo%INPUT.geo_interval!=0 or igeo<INPUT.geo_ignore) cel.read_and_used=false;
		else cel.read_and_used=true;

		cout << igeo << " " << cel.read_and_used << endl;


		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		ofs_running << "igeo=" << igeo << endl;
		cout << "igeo=" << igeo << endl;

		snapshot(ofs,cel,igeo);
	}

	ofs.close();


	return;
}

void Movie::snapshot(ofstream &ofs, const Cell &cel, const int &igeo)
{

	// get ito, ith, and itc.
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

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];
	HBs::setup_water(cel, water);



	if(INPUT.only_hydroxide==true)
	{
		int count=0;
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if(water[ia].nH!=1) continue;

			if( INPUT.nacc>0 and water[ia].naccept == INPUT.nacc)
			{
				int count=0;
				for(int iacc=0; iacc<water[ia].naccept; ++iacc)
				{
					const int iao = water[ia].acceptO[iacc];
					count += water[iao].nH;
					count += 1; // count for O
				}
				if(count!=3*INPUT.nacc) continue;
				
				int natom = (6+14)*6 + 2 + count;
				ofs << natom << endl;
				ofs << igeo << endl;

				string id="S";
				print(cel, ito, ia, ofs,id);
				for(int iadj=0; iadj<water[ia].nH; ++iadj)
				{
					const int iah = water[ia].indexH[iadj];
					print(cel, ith, iah, ofs, cel.atom[ith].id);
				}

				for(int iacc=0; iacc<water[ia].naccept; ++iacc)
				{
					const int iao = water[ia].acceptO[iacc];
					print(cel, ito, iao, ofs, cel.atom[ito].id);
					for(int iadj=0; iadj<water[iao].nH; ++iadj)
					{
						const int iah = water[iao].indexH[iadj];
						print(cel, ith, iah, ofs, cel.atom[ith].id);
					}
				}
				// print out hexane
				for(int iac=0; iac<cel.atom[itc].na; ++iac)
				{
					print(cel,itc,iac,ofs,cel.atom[itc].id);
				}
				for(int iah=191; iah<cel.atom[ith].na; ++iah)
				{
					print(cel,ith,iah,ofs,cel.atom[ith].id);
				}
			}
			else if(INPUT.nacc==-1)
			{
				int count=0;
				for(int is=0; is<INPUT.snatom; ++is)
				{
					int indexo = INPUT.satom[is]-1;
					assert(indexo>0);
					count++;
					count += water[indexo].nH;
				}

				// (6+14)*6 for hexane
				int natom = (6+14)*6+count; 
				ofs << natom << endl;
				ofs << igeo << endl;
				
				for(int is=0; is<INPUT.snatom; ++is)
				{
					int indexo = INPUT.satom[is]-1;
					print(cel, ito, indexo, ofs, cel.atom[ito].id);
					for(int iah=0; iah<water[indexo].nH; ++iah)
					{
						print(cel, ith, water[indexo].indexH[iah], ofs, cel.atom[ith].id);
					}
				}
				// print out hexane
				for(int iac=0; iac<cel.atom[itc].na; ++iac)
				{
					print(cel,itc,iac,ofs,cel.atom[itc].id);
				}
				for(int iah=191; iah<cel.atom[ith].na; ++iah)
				{
					print(cel,ith,iah,ofs,cel.atom[ith].id);
				}
			}

		}
		cout << "count=" << count << endl;
	}
	// only_hydroxide==false
	else if(INPUT.func==1)
	{
		print_special_hexane_only_ion(cel, water, ofs, ito, itc, ith, igeo);
	}
	else if(INPUT.func==2)
	{
		print_special_hexane_all_water(cel, water, ofs, itc, igeo);
	}
	else if(INPUT.func==3)
	{
		for(int is=0; is<INPUT.snatom; ++is)
		{
			int indexo = INPUT.satom[is];
			print(cel, ito, indexo, ofs, cel.atom[ito].id);
			for(int iah=0; iah<water[indexo].nH; ++iah)
			{
				print(cel, ith, water[indexo].indexH[iah], ofs, cel.atom[ith].id);
			}
		}
	}
	else if (INPUT.func == 4) 
	// print out the movie of a certain water whose accepted O atoms 
	// constitute angle less than INPUT.theta renxi 20210801
	{
		bool found = false;
		for (int ia = 0; ia < cel.atom[ito].na; ia++)
		{
			if(water[ia].naccept >= 2)
			{
				for (int iacc1 = 0; iacc1 < water[ia].naccept-1; iacc1++)
				{
					for (int iacc2 = iacc1+1; iacc2 < water[ia].naccept; iacc2++)
					{
						if (iacc1 != iacc2)
						{
							int iwater1 = water[ia].acceptO[iacc1];
							int iwater2 = water[ia].acceptO[iacc2];
							double angle = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[ia], cel.atom[ito].pos[iwater2]);
							if (angle < INPUT.theta)
							{
								found = true;
								continue;
							}
						}
					}// iacc2
					if (found == true)
					{
						continue;
					}
				}//iacc1
			}
			if (found == true)
			{
				print_within_distance(cel, ia, ofs, water, igeo, ito, ith);
			}
		}// ia
	}


	delete[] water;


	return;
}

void Movie::print(const Cell &cel, const int &it, const int &ia, ofstream &ofs, string &id)
{
	double tmpx = cel.atom[it].pos[ia].x;
	double tmpy = cel.atom[it].pos[ia].y;
	double tmpz = cel.atom[it].pos[ia].z;

	while(tmpx< 0)              tmpx+=INPUT.celldm1;
	while(tmpx>= INPUT.celldm1) tmpx-=INPUT.celldm1;
	while(tmpy< 0)              tmpy+=INPUT.celldm2;
	while(tmpy>= INPUT.celldm2) tmpy-=INPUT.celldm2;
//	if(tmpz< INPUT.z0) tmpz+=INPUT.celldm3;
	while(tmpz< 0)              tmpz+=INPUT.celldm3;
	while(tmpz>= INPUT.celldm3) tmpz-=INPUT.celldm3;

	ofs << id 
		<< " " << tmpx 
		<< " " << tmpy 
		<< " " << tmpz << endl; 
	return;
}


void Movie::print_special_hexane_only_ion(const Cell &cel, Water *water, ofstream &ofs, const int &ito, const int &itc, const int &ith, const int &igeo)
{
	int ion=-1;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(water[ia].nH==1) ion=ia;
	}
	if(ion==-1) return;
	int ntotal = (6+14)*6+2; // 2 for hydroxide 
	//ofs << INPUT.natom << endl;
	ofs << ntotal << endl;
	ofs << igeo << endl;

	int count=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			string id;
			if(it==ito and ia==ion)
			{
				//id="S";
				id="O";
				print(cel, it, ia, ofs, id);
				++count;
				for(int ih=0; ih<water[ion].nH; ++ih)
				{
					int h_index = water[ion].indexH[ih];
					print(cel, ith, h_index, ofs, cel.atom[ith].id);
					++count;
				}
			}
			else
			{
				if(it==itc)
				{
					if(ia==0 || ia==6 || ia==12 || ia==18 || ia==24 || ia==30)
					{
						id="Si";
					}
					else
					{
						id=cel.atom[it].id;
					}
					print(cel, it, ia, ofs, id);
					++count;
					for(int ih=0; ih<cel.atom[ith].na; ++ih)
					{
						double dis_ch = distance(cel.atom[itc].pos[ia], cel.atom[ith].pos[ih],
								INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						if(dis_ch < INPUT.rcut_ch)
						{
							print(cel,ith,ih,ofs,cel.atom[ith].id);
							++count;
						}
					}
				}
			}
		}//end ia
	}// end it
	cout << "count=" << count << endl;
	return;
}


void Movie::print_special_hexane_all_water(const Cell &cel, Water *water, ofstream &ofs, const int &itc, const int &igeo)
{
	ofs << INPUT.natom << endl;
	ofs << igeo << endl;

	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			if(it==itc)
			{
				string id;
				if(ia==0 || ia==6 || ia==12 || ia==18 || ia==24 || ia==30)
				{
					id="Si";
				}
				else
				{
					id=cel.atom[it].id;
				}
				print(cel, it, ia, ofs, id);
			}
			else
			{
				print(cel, it, ia, ofs, cel.atom[it].id);
			}
		}
	}
	return;
}

void Movie::print_within_distance(const Cell &cel, const int &ia, ofstream &ofs, Water* &water, const int &igeo, const int &ito, const int &ith)
{
	ofs << ">>>>>>>> This file cannot be used as VMD input. <<<<<<<" << endl;
	int nwater = 1;
	int* water_list = new int[cel.atom[ito].na];
	for (int ii = 0; ii < cel.atom[ito].na; ii++)
	{
		water_list[ii] = -1;
	}
	water_list[0] = ia;

	for (int ia1 = 0; ia1 < cel.atom[ito].na; ia1++)
	{
		if (ia1 != ia and distance(cel.atom[ito].pos[ia1], cel.atom[ito].pos[ia], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) < INPUT.rcut)
		{
			water_list[nwater] = ia1;
			nwater++;
		}
	}// ia1
	ofs << nwater*3 << endl;
	ofs << cel.snapshot_index << " " << cel.snapshot_time << " " << ia << endl;
	for (int ia2 = 0; ia2 < nwater; ia2++)
	{
		int water_index = water_list[ia2];
		Vector3<double> pos = cel.atom[ito].pos[water_index];
		double dist = special_distance(cel.atom[ito].pos[ia], pos, INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
		ofs << "O " << pos.x << " " << pos.y << " " << pos.z << endl;
		for (int iH = 0; iH < water[water_index].nH; iH++)
		{
			Vector3<double> pos1 = cel.atom[ith].pos[water[water_index].indexH[iH]];
			double dist1 = special_distance(cel.atom[ito].pos[ia], pos1, INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			ofs << "H " << pos1.x << " " << pos1.y << " " << pos.z << endl;
		}
	}
	ofs << ">>>>>>>> end <<<<<<<<" << endl;
	ofs << endl;
	return;
}

double Movie::special_distance(const Vector3<double> &pos1, Vector3<double> &pos2, const double &celldm1, const double &celldm2, const double &celldm3)
{
	double x = 0;
	double y = 0;
	double z = 0;
	if (abs(pos1.x - pos2.x) >= INPUT.celldm1)
	{
		if(abs(pos1.x - pos2.x + INPUT.celldm1) < abs(pos1.x - pos2.x) and abs(pos1.x - pos2.x + INPUT.celldm1) < abs(pos1.x - pos2.x - INPUT.celldm1))
		{
			pos2.x -= INPUT.celldm1;
			x = abs(pos1.x - pos2.x + INPUT.celldm1);
		}
		else if (abs(pos1.x - pos2.x - INPUT.celldm1) < abs(pos1.x - pos2.x) and abs(pos1.x - pos2.x - INPUT.celldm1) < abs(pos1.x - pos2.x + INPUT.celldm1))
		{
			pos2.x += INPUT.celldm1;
			x = abs(pos1.x - pos2.x - INPUT.celldm1);
		}
	}//x

	if (abs(pos1.y - pos2.y) >= INPUT.celldm2)
	{
		if(abs(pos1.y - pos2.y + INPUT.celldm2) < abs(pos1.y - pos2.y) and abs(pos1.y - pos2.y + INPUT.celldm2) < abs(pos1.y - pos2.y - INPUT.celldm2))
		{
			pos2.y -= INPUT.celldm2;
			y = abs(pos1.y - pos2.y + INPUT.celldm2);
		}
		else if(abs(pos1.y - pos2.y - INPUT.celldm2) < abs(pos1.y - pos2.y) and abs(pos1.y - pos2.y - INPUT.celldm2) < abs(pos1.y - pos2.y + INPUT.celldm2))
		{
			pos2.y += INPUT.celldm2;
			y = abs(pos1.y - pos2.y - INPUT.celldm2);
		}
	}//y

	if (abs(pos1.z - pos2.z) >= INPUT.celldm3)
	{
		if(abs(pos1.z - pos2.z + INPUT.celldm3) < abs(pos1.z - pos2.z) and abs(pos1.z - pos2.z + INPUT.celldm3) < abs(pos1.z - pos2.z - INPUT.celldm3))
		{
			pos2.z -= INPUT.celldm3;
			z = abs(pos1.z - pos2.y + INPUT.celldm2);
		}
		else if(abs(pos1.z - pos2.z - INPUT.celldm3) < abs(pos1.z - pos2.z) and abs(pos1.z - pos2.z - INPUT.celldm3) < abs(pos1.z - pos2.z + INPUT.celldm3))
		{
			pos2.z += INPUT.celldm3;
			z = abs(pos1.z - pos2.y - INPUT.celldm2);
		}
	}//z

	return sqrt(x*x+y*y+z*z);

}