#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "OH_movie.h"
#include "incremental_pdf2.h"
#include "dist2.h"
void OH_movie::Routine()
{
	int count_geometry_number = 0;

	ofstream ofs("OH.xyz");
	this->nsnapshot = 0;
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
		std::cout << "snapshot " << igeo << endl;

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
			for(int ia = 0; ia<cel.atom[ito].na; ia++)
			{
				if(water[ia].nH == 1)
				{
					if((INPUT.func == 2 and INPUT.nacc == 3 and water[ia].naccept==3) or 
						(INPUT.func == 2 and INPUT.nacc == 4 and water[ia].naccept==4) or
						(INPUT.func == 2 and INPUT.nacc == 5 and water[ia].naccept==5))
					{	
						ofs << 2+water[ia].naccept*3 << endl;
						ofs << cel.snapshot_index << " " << cel.snapshot_time << endl;
						Vector3<double>* pos = new Vector3<double>[1+water[ia].naccept*3];
						pos[0] = cel.atom[ith].pos[water[ia].indexH[0]];
						
						for (int iwater2=0; iwater2<water[ia].naccept; iwater2++)
						{
							int ia2 = water[ia].acceptO[iwater2];
							pos[1+iwater2*3] = cel.atom[ito].pos[ia2];
							pos[2+iwater2*3] = cel.atom[ith].pos[water[ia2].indexH[0]];
							pos[3+iwater2*3] = cel.atom[ith].pos[water[ia2].indexH[1]];
						}

						for (int ia_print = 0; ia_print < 1+water[ia].naccept*3; ia_print++)
						{
							incrementalPDF2::put_back(cel.atom[ito].pos[ia], pos[ia_print]);
						}

						ofs << "O " << cel.atom[ito].pos[ia].x << " " << cel.atom[ito].pos[ia].y << " " << cel.atom[ito].pos[ia].z << endl;
						ofs << "H " << pos[0].x << " " << pos[0].y << " " << pos[0].z << endl;
						for (int io = 0; io<water[ia].naccept; io++)
						{
							ofs << "O " << pos[1+io*3].x << " " << pos[1+io*3].y << " " << pos[1+io*3].z << endl;
							ofs << "H " << pos[2+io*3].x << " " << pos[2+io*3].y << " " << pos[2+io*3].z << endl;
							ofs << "H " << pos[3+io*3].x << " " << pos[3+io*3].y << " " << pos[3+io*3].z << endl;
						}
						this->nsnapshot++;
					}
				}
				if (water[ia].nH == 3)
				{
					if (INPUT.func == 1)
					{
						Vector3<double>* pos = new Vector3<double>[10];
						pos[0] = cel.atom[ito].pos[ia];
						for (int ih=0; ih<3; ih++)
						{
							pos[ih+1] = cel.atom[ith].pos[water[ia].indexH[ih]];
						}
						int ia2 = water[ia].donateO[0];
						pos[4] = cel.atom[ito].pos[ia2];
						for (int ih=0; ih<2; ih++)
						{
							pos[ih+5] = cel.atom[ith].pos[water[ia2].indexH[ih]];
						}
						int ia3 = water[ia2].donateO[0];
						pos[7] = cel.atom[ito].pos[ia3];
						for (int ih=0; ih<2; ih++)
						{
							pos[ih+8] = cel.atom[ith].pos[water[ia3].indexH[ih]];
						}
						ofs << 10 << endl;
						ofs << cel.snapshot_index << " " << cel.snapshot_time << endl;
						ofs << "O " << pos[0].x << " " << pos[0].y << " " << pos[0].z << endl;
						for (int ii=1; ii<4; ii++)
						{
							incrementalPDF2::put_back(pos[0], pos[ii]);
							ofs << "H " << pos[ii].x << " " << pos[ii].y << " " << pos[ii].z << endl;
						}
						for (int ii=0; ii<2; ii++)
						{
							
							incrementalPDF2::put_back(pos[0], pos[4+ii*3]);
							ofs << "O " << pos[4+ii*3].x << " " << pos[4+ii*3].y << " " << pos[4+ii*3].z << endl;
							for (int ih=0; ih<2; ih++)
							{
								incrementalPDF2::put_back(pos[0], pos[5+ii*3+ih]);
								ofs << "H " << pos[5+ii*3+ih].x << " " << pos[5+ii*3+ih].y << " " << pos[5+ii*3+ih].z << endl;
							}
						}
						this->nsnapshot++;
					}
					else if (INPUT.func == 2 and water[ia].ndonate == 3 and water[water[ia].donateO[0]].ndonate >= 1)
					{
						Vector3<double>* pos = new Vector3<double>[16];
						pos[0] = cel.atom[ito].pos[ia];
						ofs << 16 << endl;
						ofs << "O " << pos[0].x << " " << pos[0].y << " " << pos[0].z << endl;
						for (int ih=0; ih<3; ih++)
						{
							pos[ih+1] = cel.atom[ith].pos[water[ia].indexH[ih]];
							incrementalPDF2::put_back(pos[0], pos[ih+1]);
							ofs << "H " << pos[ih+1].x << " " << pos[ih+1].y << " " << pos[ih+1].z << endl;
						}
						for (int idon=0; idon<3; idon++)
						{
							pos[4+idon*3] = cel.atom[ito].pos[water[ia].donateO[idon]];
							pos[4+idon*3+1] = cel.atom[ith].pos[water[water[ia].donateO[idon]].indexH[0]];
							pos[4+idon*3+2] = cel.atom[ith].pos[water[water[ia].donateO[idon]].indexH[1]];
							incrementalPDF2::put_back(pos[0], pos[4+idon*3]);
							ofs << "O " << pos[4+idon*3].x << " " << pos[4+idon*3].y << " " << pos[4+idon*3].z << endl;
							incrementalPDF2::put_back(pos[0], pos[5+idon*3]);
							ofs << "H " << pos[5+idon*3].x << " " << pos[5+idon*3].y << " " << pos[5+idon*3].z << endl;
							incrementalPDF2::put_back(pos[0], pos[6+idon*3]);
							ofs << "H " << pos[6+idon*3].x << " " << pos[6+idon*3].y << " " << pos[6+idon*3].z << endl;
						}
						int second_neighbor_index = water[water[ia].donateO[0]].donateO[0];
						pos[13] = cel.atom[ito].pos[second_neighbor_index];
						incrementalPDF2::put_back(pos[0], pos[13]);
						ofs << "O " << pos[13].x << " " << pos[13].y << " " << pos[13].z << endl;
						pos[14] = cel.atom[ith].pos[water[second_neighbor_index].indexH[0]];
						incrementalPDF2::put_back(pos[0], pos[14]);
						ofs << "H " << pos[14].x << " " << pos[14].y << " " << pos[14].z << endl;
						pos[15] = cel.atom[ith].pos[water[second_neighbor_index].indexH[1]];
						incrementalPDF2::put_back(pos[0], pos[15]);
						ofs << "H " << pos[15].x << " " << pos[15].y << " " << pos[15].z << endl;
						this->nsnapshot++;
					}
					else if (INPUT.func == 4 and water[ia].naccept==1 and water[water[ia].donateO[0]].naccept == 2) // plot the presolvation structure of hydronium
					{
						Vector3<double>* pos = new Vector3<double>[13];
						pos[0] = cel.atom[ito].pos[ia];
						ofs << 13 << endl;
						ofs << "O " << pos[0].x << " " << pos[0].y << " " << pos[0].z << endl;
						for (int ih=0; ih<3; ih++)
						{
							pos[ih+1] = cel.atom[ith].pos[water[ia].indexH[ih]];
							incrementalPDF2::put_back(pos[0], pos[ih+1]);
							ofs << "H " << pos[ih+1].x << " " << pos[ih+1].y << " " << pos[ih+1].z << endl;
						}
						int O1 = water[ia].acceptO[0];
						pos[4] = cel.atom[ito].pos[O1];
						incrementalPDF2::put_back(pos[0], pos[4]);
						ofs << "O " << pos[4].x << " " << pos[4].y << " " << pos[4].z << endl;
						for (int ih=0; ih<2; ih++)
						{
							pos[ih+5] = cel.atom[ith].pos[water[O1].indexH[ih]];
							incrementalPDF2::put_back(pos[4], pos[ih+5]);
							ofs << "H " << pos[ih+5].x << " " << pos[ih+5].y << " " << pos[ih+5].z << endl;
						}
						int O2 = water[ia].donateO[0];
						pos[7] = cel.atom[ito].pos[O2];
						incrementalPDF2::put_back(pos[0], pos[7]);
						ofs << "O " << pos[7].x << " " << pos[7].y << " " << pos[7].z << endl;
						for (int ih=0; ih<2; ih++)
						{
							pos[ih+8] = cel.atom[ith].pos[water[O2].indexH[ih]];
							incrementalPDF2::put_back(pos[7], pos[ih+8]);
							ofs << "H " << pos[ih+8].x << " " << pos[ih+8].y << " " << pos[ih+8].z << endl;
						}
						int O3 = -1;
						for (int iacc=0; iacc<2; iacc++)
						{
							if (water[O2].acceptO[iacc] != ia) O3 = water[O2].acceptO[iacc];
						}
						pos[10] = cel.atom[ito].pos[O3];
						incrementalPDF2::put_back(pos[0], pos[10]);
						ofs << "O " << pos[10].x << " " << pos[10].y << " " << pos[10].z << endl;
						for (int ih=0; ih<2; ih++)
						{
							pos[ih+11] = cel.atom[ith].pos[water[O3].indexH[ih]];
							incrementalPDF2::put_back(pos[10], pos[ih+11]);
							ofs << "H " << pos[ih+11].x << " " << pos[ih+11].y << " " << pos[ih+11].z << endl;
						}
						this->nsnapshot++;
					} // func == 4
				}
			}
			if (INPUT.func == 3)
			{
				assert(INPUT.n_recorded_water > 0);
				assert(INPUT.Oindex[0] >= 0);
				if (INPUT.system == "hydroxide") ofs << 3*INPUT.n_recorded_water-1 << endl;
				else if (INPUT.system == "hydronium") ofs << 3*INPUT.n_recorded_water+1 << endl;
				else if (INPUT.system == "water") ofs << 3*INPUT.n_recorded_water << endl;
				ofs << cel.snapshot_index << " " << cel.snapshot_time << endl;
				Vector3<double> initial_pos; 
				for (int io=0; io<INPUT.n_recorded_water; io++)
				{
					int water_index = INPUT.Oindex[io]-1;
					assert(water_index >= 0);
					if (io == 0)
					{
						initial_pos = cel.atom[ito].pos[water_index];
						ofs << "O " << initial_pos.x << " " << initial_pos.y << " " << initial_pos.z << endl;
					}
					if (io != 0)
					{
						Vector3<double> instant_pos = cel.atom[ito].pos[water_index];
						Dist2::putback_cell(initial_pos, instant_pos);
						ofs << "O " << instant_pos.x << " " << instant_pos.y << " " << instant_pos.z << endl;
					}
					for (int ih=0; ih<water[water_index].nH; ih++)
					{
						Vector3<double> instant_pos = cel.atom[ith].pos[water[water_index].indexH[ih]];
						Dist2::putback_cell(initial_pos, instant_pos);
						ofs << "H " << instant_pos.x << " " << instant_pos.y << " " << instant_pos.z << endl;
					}
				}
			}
		}
		if (this->nsnapshot >= INPUT.ntry)
		{
			break;
		}
		delete[] water;
		cel.clean();
	}
	ofs.close();
}


