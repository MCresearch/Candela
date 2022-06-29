#include "cellFile.h"
#include "input.h"
#include "tetra_order.h"
#include "math.h"
#include "HBs.h"

void TOP::Routine()
{
	TITLE("TOP","Routine");
	
	cout << "Compute the Tetrahedral Order Parameter for Water" << endl;

	count_adj=0;
	count_adj2=0;
	tetra_order_parameter=0.0;

	this->dr = INPUT.dr;
	assert(dr>0.0);
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);
	this->nmesh = int(rcut / dr) +  1;
	
	this->dis_tetra = new double[nmesh]();
	this->dis_count = new int[nmesh]();
	

	// TOP (tetrahedral order parameter) in terms of water density
	this->dr_den=0.005; // g/cm^3 is the unit
	this->dcut=1.20; // g/cm^3 is the unit
	this->nden = int(dcut/dr_den)+1;	
	this->den_tetra = new double[nden]();
	this->den_count = new int[nden]();

	//**************************
	// BEGIN CALCULATING DATA
	//**************************
	int count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) continue;

		// cel : input geometry file
		CellFile cel;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;

		// calculate the density of liquid water (only for 64 water molecules)
		this->density = 64*18*1.6605/cel.volume;
		cout << "igeo=" << igeo << " " << density << " g/cm^3" << endl;

		cal_tetra_order_para(cel);
	}	

	if(count_adj>0)
	{
		tetra_order_parameter=tetra_order_parameter/(double)count_adj;
		ave_adj=ave_adj/(double)count_adj2;
		cout << "tetra_order_parameter = " << tetra_order_parameter << endl;
		cout << "average adjacent O atoms = " << setprecision(6) << ave_adj << endl;
	} 

	if(INPUT.func_b==2)
	{
		ofstream ofs("dis_tetra.dat");
		for(int ir=0; ir<nmesh; ++ir)
		{
			if(dis_count[ir]>0)
			{
				ofs << ir * dr << " " << dis_tetra[ir]/dis_count[ir] << endl;
			}
		}
	}

	delete[] dis_tetra;
	delete[] dis_count;

	// print out the TOP in terms of water density
	ofstream ofs_d("den_tetra.dat");
	for(int i=0; i<nden; ++i)
	{
		if(den_count[i]>0)
		{
			ofs_d << i * this->dr_den << " " << den_tetra[i]/den_count[i] << endl;
		}
	}
	ofs_d.close();

	return;
}


void TOP::cal_tetra_order_para(const Cell &cel)
{
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

	// search for O-O pair
	int ind_o[12];
	double dis_o[12];
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
// no need for this function, mohan added 2017-03-05
//		if( cel.atom[ito].pos[ia].z > INPUT.z0 and cel.atom[ito].pos[ia].z < INPUT.z1 )
		{
			// step 1: initialize 
			int iadj=0;
			for(int i=0; i<12; ++i){ ind_o[i]=-1;dis_o[i]=0.0;}

			// step 2: distance
			for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
			{
				if(ia==ia2) continue;
				double dis = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				if(dis < INPUT.bdf_rcut)
				{
					// test
					if(iadj>=4)
					{
						for(int iii=0; iii<4; ++iii)
						{
							if(dis<dis_o[iii])
							{
								ind_o[iii]=ia2;
								dis_o[iii]=dis;
								break;
							}		
						}
					}
					else
					{
						ind_o[iadj]=ia2;
						dis_o[iadj]=dis;
					}
					++iadj;
				}
			}

			if(iadj>=4)
			{
				double value=0.0;	
				for(int i=0; i<3; ++i)
				{
					for(int j=i+1; j<4; ++j)
					{
						double angle0 = HBs::angle(cel, cel.atom[ito].pos[ind_o[i]], cel.atom[ito].pos[ia], cel.atom[ito].pos[ind_o[j]]);
						value += pow((cos(angle0/180*3.1415926535)+1.0/3.0), 2);
						//					cout << ind_o[i] << " " << ia << " " << ind_o[j] << endl;
						//					cout << "v=" << value << " angle0=" << angle0 << endl;
					}
				}
				value = value*3.0/8.0;

				// this is for surface water regarding the instatanous liquid interface
				if(INPUT.func_b==2)
				{
					int index_jj = cel.atom[ito].pos_ili[ia]/this->dr;
					if(index_jj>=0 and index_jj<nmesh)
					{
						dis_tetra[index_jj] += 1-value;
						dis_count[index_jj] += 1;
					}
				}

				tetra_order_parameter += 1-value;
				
				// ion density
				int index_den=this->density/this->dr_den;
				this->den_tetra[index_den] += 1-value;
				this->den_count[index_den] += 1;

				ofs_running << "tetra " << cel.snapshot_time << " " << 1-value << endl;
				++count_adj;
			}

			ave_adj+=iadj;
			++count_adj2;

		}
	}

	return;
}


