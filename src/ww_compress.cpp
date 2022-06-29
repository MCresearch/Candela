#include "cellFile.h"
#include "input.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"
#include "ww_compress.h"

WW_Compress::WW_Compress()
{
}

WW_Compress::~WW_Compress(){}

void WW_Compress::Routine()
{
	TITLE("Waterwire","Routine");
	
	cout << "Compute the compression of waterwires" << endl;

	this->dr = INPUT.dr;
	assert(dr>0.0);

	// radius cutoff in real space, usually choose a/2,
	// where a is the lattice constant.
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);

	// number of radial mesh.
	this->nmesh = int(rcut / dr) +  1;

	this->ww = new double[nmesh]();
	this->ww_all = new double[nmesh]();
	this->ww1 = new double[nmesh]();
	this->ww2 = new double[nmesh]();

	this->ww_h2o = new double[nmesh]();
	this->ww_h2o_all = new double[nmesh]();
	this->ww1_h2o = new double[nmesh]();
	this->ww2_h2o = new double[nmesh]();




	// BEGIN CALCULATING DATA
	this->count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

//		if(cel.snapshot_time > 1.0) continue; //mohan test

		if(cel.read_and_used==false) continue;
		cout << "igeo=" << igeo << endl;
		count_compression(cel, igeo);
		++count_geometry_number;
	}	



	// print out data
	ofstream ofs_result("waterwire.dat");

	double sum=0.0;
	double sum1=0.0;
	double sum2=0.0;
	for(int i=0; i<nmesh; ++i) sum += ww[i]*dr;
	for(int i=0; i<nmesh; ++i) sum1 += ww1[i]*dr;
	for(int i=0; i<nmesh; ++i) sum2 += ww2[i]*dr;

	double sumw=0.0;
	double sum1w=0.0;
	double sum2w=0.0;
	double sum_all=0.0;
	for(int i=0; i<nmesh; ++i) sumw += ww_h2o[i]*dr;
	for(int i=0; i<nmesh; ++i) sum1w += ww1_h2o[i]*dr;
	for(int i=0; i<nmesh; ++i) sum2w += ww2_h2o[i]*dr;
	for(int i=0; i<nmesh; ++i) sum_all += ww_h2o_all[i]*dr;

//	if(sum>0.0 and sum1>0.0 and sum2>0.0)
	{
		ofs_result << "r sum_ion sum_dis1 sum_dis2 sum_water sum_dis1 sum_dis2 free_water" << endl;
		double ss1=0;
		double ss2=0;
		for(int i=0; i<nmesh; ++i)
		{
//			double free_energy = -std::log(ww_h2o[i]*INPUT.factor); 
//			ofs_result << i*dr << " " << ww[i]/sum << " " << ww1[i]/sum1 << " " << ww2[i]/sum2 << " "
//				<< ww_h2o[i]/sumw << " " << ww1_h2o[i]/sum1w << " " << ww2_h2o[i]/sum2w << " " << free_energy << endl;

//			ofs_result << i*dr << " " << ww_h2o[i]/(double)count_geometry_number  << " " << ww_h2o[i]/sum
//					<< " " << ww_h2o_all[i]/(double)count_geometry_number << " " << ww_h2o_all[i]/sum_all << endl; 

			ofs_result << i*dr << " " << ww_h2o[i]/ww_h2o_all[i] << endl;

			ss1 += ww_h2o[i];
			ss2 += ww_h2o_all[i];

//			ss1 += ww[i];
//			ss2 += ww_all[i];
		}
		if(ss2>0)
		{
			ofs_running << "count_geometry_number" << count_geometry_number << endl;
			ofs_running << "sum all HB-connected h2o waterwire is " << ss1/(double)count_geometry_number/64.0 
			<< " ratiio is " << ss1/ss2 << endl;
			ofs_running << "sum all possible h2o waterwire is " << ss2/(double)count_geometry_number/64.0 << endl;
		}
	}	
	ofs_result.close();


	delete[] ww;
	delete[] ww1;
	delete[] ww2;

	delete[] ww_h2o;
	delete[] ww1_h2o;
	delete[] ww2_h2o;

	delete[] ww_all;
	delete[] ww_h2o_all;

	return;
}


void WW_Compress::count_compression(const Cell &cel, const int &igeo)
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


	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

	HBs::setup_water(cel, water);
	ofs_running << "snapshot " << cel.snapshot_index << " nions = " << Water::nions << endl;

	neighbour1(cel, water, ito, ith);

	delete[] water;
}

void WW_Compress::neighbour1(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		//cout << "nH=" << water[ia].nH << endl;
//		if(water[ia].nH==2)
		{
			// all possible water wires
			for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
			{
//				if(water[ia2].nH!=2) continue;
				if(ia==ia2) continue;

/*
				bool firstHB=false;
				for(int iacc=0; iacc<water[ia].ndonate; ++iacc)
				{
					int io_adj = water[ia].donateO[iacc];
					if(ia2==io_adj) firstHB=true;
				}
				for(int iacc=0; iacc<water[ia].naccept; ++iacc)
				{
					int io_adj = water[ia].acceptO[iacc];
					if(ia2==io_adj) firstHB=true;
				}
				if(firstHB==false)continue;
*/
				double dis1 = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				if(dis1 >= INPUT.rcut) continue;
//				if(dis1 > 2.75) continue;
				for(int ia3=0; ia3<cel.atom[ito].na; ++ia3)
				{
//					if(water[ia3].nH!=2) continue;
					if(ia3==ia2) continue;
					if(ia3==ia) continue;
					double dis2 = distance(cel.atom[ito].pos[ia2], cel.atom[ito].pos[ia3], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
//					if(dis1 <= dis2)
					{
						double dis12 = dis1 + dis2;
						if(dis12 < INPUT.rcut)
						{
							int index = dis12/dr; ww_h2o_all[index]++;
							//cout << ww_h2o_all[index] << endl;
						}
					}
				}
			}

			// donate + donate waterwire
			for(int iacc=0; iacc<water[ia].ndonate; ++iacc)
			{
				double dis1 = water[ia].donate_disO[iacc];
				int io_adj = water[ia].donateO[iacc];
//				if(water[io_adj].nH!=2) continue;
				
				for(int iacc2=0; iacc2<water[io_adj].ndonate; ++iacc2)
				{
					if( water[io_adj].donateO[iacc2]==ia ) continue;
					if( water[io_adj].donateO[iacc2]==io_adj ) continue;

					int io_adj2 = water[io_adj].donateO[iacc2];
//					if(water[io_adj2].nH!=2) continue;

					double dis2 = water[io_adj].donate_disO[iacc2];
					double dis12 = dis1 + dis2;
					
					if(dis12 < INPUT.rcut)
					{
						int index = dis12/dr; ww_h2o[index]++;
						index = dis1/dr; ww1_h2o[index]++;
						index = dis2/dr; ww2_h2o[index]++;
					}
				}
			}

		}

		if(INPUT.system=="hydroxide" and water[ia].nH==1 and Water::nions==1)
		{
			for(int iacc=0; iacc<water[ia].naccept; ++iacc)
			{
				double dis1 = water[ia].accept_disO[iacc];
				int io_adj = water[ia].acceptO[iacc];
				
				for(int iacc2=0; iacc2<water[io_adj].naccept; ++iacc2)
				{
					double dis2 = water[io_adj].accept_disO[iacc2];
					double dis12 = dis1 + dis2;
					
					if(dis12 < INPUT.rcut)
					{
						//if(water[ia].naccept==3)
						{
							int index = dis12/dr; ww[index]++;
							index = dis1/dr; ww1[index]++;
							index = dis2/dr; ww2[index]++;
						}
					}
				}
			}
		} 
		else if(INPUT.system=="hydronium" and water[ia].nH==3 and Water::nions==1)
		{
			// possible waterwire for hydronium
			for(int idon=0; idon<water[ia].ndonate; ++idon)
			{
				double dis1 = water[ia].donate_disO[idon];
				int io_adj = water[ia].donateO[idon];
				
				for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
				{
					double dis2 = distance(cel.atom[ito].pos[ia2], cel.atom[ito].pos[io_adj], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					double dis12 = dis1 + dis2;
	
					if(dis12 < INPUT.rcut)
					{
						int index = dis12/dr; ww_all[index]++;
					}
				}
			}

			// possible HB-connected waterwire for hydronium
			for(int idon=0; idon<water[ia].ndonate; ++idon)
			{
				double dis1 = water[ia].donate_disO[idon];
				int io_adj = water[ia].donateO[idon];
				
				for(int idon2=0; idon2<water[io_adj].ndonate; ++idon2)
				{
					double dis2 = water[io_adj].donate_disO[idon2];
					double dis12 = dis1 + dis2;
	
					if(dis12 < INPUT.rcut)
					{
						int index = dis12/dr; ww[index]++;
						index = dis1/dr; ww1[index]++;
						index = dis2/dr; ww2[index]++;
					}
				}
			}
		}
	}
	return;
}
