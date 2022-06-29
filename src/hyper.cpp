#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "hyper.h"
#include "math.h"
#include "water.h"
#include "pdf_added.h" // to use compute_delta
#include "mj.h"

Hyper::Hyper()
{
	ave_donate=0.0;
}

Hyper::~Hyper(){}

void Hyper::Routine()
{
	cout << "compute order-parameters for tetra and hyper coordinates"  << endl;

	this->dr = INPUT.dr;
	assert(dr>0.0);
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);
	// number of radial mesh.
	this->nmesh = int(rcut / dr) +  1;
	
	double* gr = new double[nmesh]();


	// mohan added 2016-12-03
	ProtonTransfer PT;
	PT.setup_PT();

	
	this->count_geometry_number=0;
	for(this->igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(PT.use_pt==true)
		{
			string type=PT.which_type_pt(cel.snapshot_index);
			ofs_running << "snapshot " << cel.snapshot_index << " " << type << endl;
			cout << "snapshot " << cel.snapshot_index << " " << type << endl;
			if(INPUT.func_b==2) // for rattling (rest)
			{
				if(type!="rattling") continue;
			}
			else if(INPUT.func_b==3) // for all PT
			{
				if(type=="rattling") continue;
			}
			else if(INPUT.func_b==4) // for multiple jumps
			{
				if(type!="double" and type!="triple" and type!="quadralple" and type!="five") continue;
			}
		}

		if(cel.read_and_used==false) continue;
		cout << "igeo=" << igeo << endl;

		hypercoordinated(cel, gr);
	}

	assert(count_geometry_number>0);
	this->ave_donate = this->ave_donate/(double)count_geometry_number;
	cout << "average donate HB number is " << ave_donate << endl;
	ofs_running << "average donate HB number is " << ave_donate << endl;

	ofstream ofs(INPUT.geo_out.c_str());
	// we output the pair correlation function and static structrue factor
	double sum=0.0;
	for(int i=0; i<nmesh-1; ++i)
	{
		gr[i] = gr[i]/(double)count_geometry_number/dr;
		sum += gr[i]*dr;
		ofs << (i+0.5)*dr << " " << gr[i] << endl;
	}
	cout << "count events: " << count_geometry_number << endl;
	cout << "sum is " << sum << endl;
	ofs.close();

	delete[] gr;
	
	return;
}


void Hyper::hypercoordinated(const Cell &cel, double * gr)
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
	Water::nions=0;

	HBs::setup_water(cel, water);

	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(INPUT.system=="hydroxide" and water[ia].nH==1)
		{
			assert(INPUT.nacc>0);

			// case 3 acccepted and 1 donte
			if(water[ia].naccept==INPUT.nacc and water[ia].ndonate==INPUT.ndon)
			{
				one_donate(cel, water, ito, ia, gr);
			}
			// case 4 accepted
			else if(water[ia].naccept==INPUT.nacc and INPUT.ndon==-1)
			{
				assert(INPUT.nacc==4);
				if(INPUT.delta!=0.0)
				{
					double delta = PDF_ADDED::compute_delta(cel,water,ito,ia);
					cout << "delta " << delta << endl;
					if(INPUT.delta<0.0)
					{
						if(delta>=abs(INPUT.delta)) continue;
					}
					else if(INPUT.delta>0.0)
					{
						if(delta<=abs(INPUT.delta)) continue;
					}
				}
				if(INPUT.func_c==2 and water[ia].ndonate==0)
				{
					continue;
				} 
				if(INPUT.func_c==3 and water[ia].ndonate != 0)
				{
					continue;
				} // renxi added 20201031
				// func_c = 2: calculate planarity of Hydroxide with donation != 0
				// func_c = 3: calculate planarity of Hydroxide with donation =0
				four_accept(cel, water, ito, ia, gr);
			}
			else if(INPUT.nacc==30)
			{
				if(INPUT.delta!=0.0)
				{
					double delta = PDF_ADDED::compute_delta(cel,water,ito,ia);
					cout << "delta " << delta << endl;
					if(INPUT.delta<0.0)
					{
						if(delta>=abs(INPUT.delta)) continue;
					}
					else if(INPUT.delta>0.0)
					{
						if(delta<=abs(INPUT.delta)) continue;
					}
				} 
				if(water[ia].naccept==3 and water[ia].ndonate==1)
				{
					one_donate(cel, water, ito, ia, gr);
				}
				if(water[ia].naccept==4)
				{
					four_accept(cel, water, ito, ia, gr);
				}
			}// nacc 30
		}
	}

	delete[] water;
	
	return;
}

void Hyper::dvec(const Cell &cel, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &da)
{
	double celldm1_half = cel.a1.norm()/2.0;
	double celldm2_half = cel.a2.norm()/2.0;
	double celldm3_half = cel.a3.norm()/2.0;

	double dx, dy, dz;

	dx = a1.x - a2.x;
	while( dx >  celldm1_half ) dx -= cel.a1.norm(); 
	while( dx < -celldm1_half ) dx += cel.a1.norm(); 

	dy = a1.y - a2.y;
	while( dy >  celldm2_half ) dy -= cel.a2.norm(); 
	while( dy < -celldm2_half ) dy += cel.a2.norm(); 

	dz = a1.z - a2.z;
	while( dz >  celldm3_half ) dz -= cel.a3.norm(); 
	while( dz < -celldm3_half ) dz += cel.a3.norm(); 

	da.x = dx;
	da.y = dy;
	da.z = dz;
	

// for tests only
//	ofs_running << endl;
//	ofs_running << "dx dy dz are " << dx << " " << dy << " " << dz << endl;
//	ofs_running << "disance is " << sqrt(dx*dx+dy*dy+dz*dz) << " angstroms" << endl;

	return;
} 


void Hyper::one_donate(const Cell &cel, const Water* water, const int &ito, const int &ia, double* gr)
{
	Vector3<double> da1, da2, da3, n, d;
	assert(INPUT.nacc==3 or INPUT.nacc==30);
	const int io1 = water[ia].donateO[0];
	const int io2 = water[ia].acceptO[0];
	const int io3 = water[ia].acceptO[1];
	const int io4 = water[ia].acceptO[2];

	Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io2], da1);
	Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io3], da2);
	Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io4], da3);

	n = da1^da2;
	n = n/n.norm();
	d = n*da3;
	double dis = d.norm(); 
	const int ind_grid = (double)dis/dr; 
	if(ind_grid<nmesh)
	{
		gr[ind_grid] = gr[ind_grid]+1.0;
		ave_donate+=1.0;
		++count_geometry_number;
		//	cout << count_geometry_number << endl;
	}



	return;
}

void Hyper::four_accept(const Cell &cel, const Water *water, const int &ito, const int &ia, double* gr)
{
	Vector3<double> da1, da2, da3, n, d;
	double average=0.0;

	ofs_running << "atom list: " << ia << " " << water[ia].acceptO[0]
		<< " " << water[ia].acceptO[1] << " " << water[ia].acceptO[2]
		<< " " << water[ia].acceptO[3] << endl;

	ofs_running << cel.snapshot_index << " " << cel.snapshot_time; 


	for(int i1=0; i1<water[ia].naccept; ++i1)
	{
		int i2=(i1+1)%4;
		int i3=(i1+2)%4;
		int i4=(i1+3)%4;

		const int io1 = water[ia].acceptO[i1];
		const int io2 = water[ia].acceptO[i2];
		const int io3 = water[ia].acceptO[i3];
		const int io4 = water[ia].acceptO[i4];

		Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io2], da1);
		Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io3], da2);
		Hyper::dvec(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io4], da3);

		n = da1^da2;
		n = n/n.norm();
		d = n*da3;
		double dis = d.norm(); 

		ofs_running << " " << dis;

		if(INPUT.func==1) // record for all X
		{
			const int ind_grid = (double)dis/dr; 
			if(ind_grid<nmesh)
			{
				gr[ind_grid] = gr[ind_grid]+1.0;
				ave_donate+=water[ia].ndonate;
				++count_geometry_number;
				//	cout << count_geometry_number << endl;
			}
		}
		else if(INPUT.func==2) // only record for the average X
		{
			average += dis;
		}
	}

	ofs_running << endl;

	if(INPUT.func==1)
	{
		// do nothing
	}
	else if(INPUT.func==2)
	{
		average/=4.0;
		ofs_running << cel.snapshot_index << " " << igeo << " average " << average << endl; 
		const int ind_grid = average/dr;
		if(ind_grid<nmesh)
		{
			gr[ind_grid] = gr[ind_grid]+1.0;
			ave_donate+=water[ia].ndonate;
			++count_geometry_number;
		}
	}


//	ofs_running << igeo << " order_x " << average << endl;
	return;
}
