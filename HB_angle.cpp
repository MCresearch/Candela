#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "HB_angle.h"

//HB_angle::HB_angle(){}

//HB_angle::~HB_angle(){}

void HB_angle::Routine()
{
	int count_geometry_number = 0;
	double* angles;
	double* dis_OH2_list;
	if (INPUT.func == 1 or INPUT.func == 2)
	{
		angles = new double[10]();
		dis_OH2_list = new double[10]();
	}
	else if (INPUT.func == 3)
	{
		angles = new double[int(INPUT.theta/INPUT.dtheta)]();
		for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
		{
			angles[ia] = 0;
		}
	}
	double dist_OH = -1;
	double dist_HH = -1;
	double dist_OH2 = -1;
	double average_angle = -1;
	double average_dis = -1;
	ofstream ofs("angles.txt");
	ofstream ofs_dis("dis.txt");
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

		if(cel.read_and_used==false) 
		{
			cel.clean();//qianrui add in 2020-1-7
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;

		int ito = -1;
		int ith = -1;
		for(int it=0; it<INPUT.ntype; it++)
		{
			if(cel.atom[it].id=="O"){ito = it;}
			else if(cel.atom[it].id=="H" or cel.atom[it].id=="D"){ith = it;}
		}
		assert(ito>=0 and ith>=0);
		Water *water = new Water[cel.atom[ito].na];
		Water::nions = 0;
		HBs::setup_water(cel, water);
		if(Water::nions == 1)
		{
			if(INPUT.system == "hydroxide")
			{
				for(int ia = 0; ia<cel.atom[ito].na; ia++)
				{
					if(water[ia].nH == 1)
					{
						if (INPUT.func == 1 or INPUT.func == 2)
						{
							for(int iH=0; iH<10; iH++)
							{
								angles[iH] = 0;
								dis_OH2_list[iH] = 0;
							}
							dist_OH = distance(cel.atom[ito].pos[ia], cel.atom[ith].pos[water[ia].indexH[0]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							dist_OH2 = -1;
							dist_HH = -1;
							average_angle = 0;
							average_dis = 0;
							for(int iH = 0; iH < water[ia].naccept; iH++)
							{
								if (INPUT.func == 1)
								{
									dist_OH2 = distance(cel.atom[ito].pos[ia], cel.atom[ith].pos[water[ia].acceptH[iH]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
									dist_HH = distance(cel.atom[ith].pos[water[ia].indexH[0]], cel.atom[ith].pos[water[ia].acceptH[iH]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
								}
								else if (INPUT.func == 2)
								{
									dist_OH2 = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[water[ia].acceptO[iH]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
									dist_HH = distance(cel.atom[ith].pos[water[ia].indexH[0]], cel.atom[ito].pos[water[ia].acceptO[iH]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
								}
								angles[iH] = acos((pow(dist_OH2, 2) + pow(dist_OH, 2) - pow(dist_HH, 2))/2/dist_OH/dist_OH2)/3.1415926535897*180;
								average_angle += angles[iH];
								dis_OH2_list[iH] = dist_OH2;
								average_dis += dist_OH2;
							}
							average_angle /= water[ia].naccept;
							average_dis /= water[ia].naccept;
							ofs << cel.snapshot_index << " " << cel.snapshot_time << " " << average_angle;
							ofs_dis << cel.snapshot_index << " " << cel.snapshot_time << " " << average_dis;
							for(int iH=0; iH<10; iH++)
							{
								ofs << " " << angles[iH];
								ofs_dis << " " << dis_OH2_list[iH];
							}
							ofs << endl;
							ofs_dis << endl;
						}
						else if (INPUT.func == 3)
						{
							for (int ia2=0; ia2<cel.atom[ito].na; ia2++)
							{
								if (ia2 == ia) continue;
								double dis = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
								if (dis > INPUT.rcut) continue;
								int Hindex = water[ia].indexH[0];
								double angle = HBs::angle(cel, cel.atom[ith].pos[Hindex], cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2]);
								if (angle < INPUT.theta)
								{
									int which = int(angle/INPUT.dtheta);
									angles[which]++;
								}
							}
						}
					}
				}
			}
			else if (INPUT.system == "hydronium")
			{
				for(int ia=0; ia<cel.atom[ito].na; ia++)
				{
					if(water[ia].nH == 3)// locate hydronium
					{
						for(int iH=0; iH<10; iH++)
						{
							angles[iH] = 0;
							dis_OH2_list[iH] = 0;
						}
						dis_OH2_list[0] = water[ia].naccept;
						angles[0] = water[ia].naccept; // ATTENTION! Here angles restore length of O-H covalent bond in hydronium.
						for(int ih=0; ih<water[ia].ndonate; ih++)
						{
							int indexH = water[ia].donateH[ih];
							int indexO = water[ia].donateO[ih];
							angles[ih+1] = distance(cel.atom[ito].pos[ia], cel.atom[ith].pos[indexH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							dis_OH2_list[ih+1] = distance(cel.atom[ito].pos[indexO], cel.atom[ith].pos[indexH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						}
						ofs_dis << cel.snapshot_index << " " << cel.snapshot_time << " ";
						ofs << cel.snapshot_index << " " << cel.snapshot_time << " ";
						for(int ih=0; ih<10; ih++)
						{
							ofs << angles[ih] << " ";
							ofs_dis << dis_OH2_list[ih] << " ";

						}
						ofs << endl;
						ofs_dis << endl;
					}
				}
			}
		}
	}
	if (INPUT.func == 3)
	{
		double sum = 0;
		for (int ia=0; ia < int(INPUT.theta/INPUT.dtheta); ia++)
		{
			sum += angles[ia]*INPUT.dtheta;
		}
		for (int ia=0; ia < int(INPUT.theta/INPUT.dtheta); ia++)
		{
			angles[ia] /= sum;
			ofs << (ia + 0.5)*INPUT.dtheta << " " << angles[ia] << endl;
		}
	}
	ofs.close();
	ofs_dis.close();
}
