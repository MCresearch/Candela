#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "trajadj.h"
#include "math.h"
#include "water.h"
#include "pdf.h" // to use compute_delta

Trajadj::Trajadj(){}
Trajadj::~Trajadj(){}

void Trajadj::Routine()
{
	cout << "compute the traj-adj 2D figure"  << endl;

	assert(INPUT.ntzone>0);

	this->dr = INPUT.dr;
	assert(dr>0.0);
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);
	this->nmesh = int(rcut / dr)+1;

	cout << "ntzone = " << INPUT.ntzone << endl;
	cout << "nmesh = " << nmesh << endl;

	double** trmap = new double*[INPUT.ntzone];
	for(int i=0;i<INPUT.ntzone;++i)
	{
		trmap[i] = new double[nmesh]();
	}

	cout << "the trajectory is: " << endl;
	cout << INPUT.geo_1 << " " << INPUT.geo_2 << endl;
	
	int ngeom = INPUT.geo_2-INPUT.geo_1+1;
	cout << "the number of geometries is: " << ngeom << endl;
	int ngeom_each = ngeom/INPUT.ntzone;
	cout << "the number of geometry in each sector: " << ngeom_each << endl;
	if(ngeom%INPUT.ntzone!=0)
	{
		cout << "ngeom%INPUT.ntzone==" << ngeom%INPUT.ntzone << endl;
		exit(0);
	}

	assert(INPUT.geo_1>0);

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		cout << "igeo=" << igeo << endl;
		this->itzone = (igeo-1)/ngeom_each;
//		cout << "in time zone: " << itzone << endl;

		compute(cel, trmap, igeo);
	}

	
	ofstream ofs(INPUT.geo_out.c_str());
	for(int i=0; i<INPUT.ntzone; ++i)
	{
		for(int j=0; j<nmesh-1; ++j)
		{
			ofs << trmap[i][j] << " ";
		}
		ofs << endl;
	}
	ofs.close();

	// clean
	for(int i=0; i<INPUT.ntzone; ++i)
	{
		delete[] trmap[i];
	}
	delete[] trmap;

	return;
}

void Trajadj::compute(const Cell &cel, double** trmap, const int &igeo)
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

	// search for OH
	int ion=-1;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
		{
			double ho_dis = distance(cel.atom[ito].pos[ia],cel.atom[ith].pos[ia2],norm1,norm2,norm3);
			if(ho_dis<INPUT.rcut_oh)
			{
				int ind=water[ia].nH;
				water[ia].indexH[ind]=ia2;
				water[ia].disH[ind]=ho_dis;
				water[ia].nH++;	
			}
		}
		if(water[ia].nH==1) ion=ia;
	}


	for(int io=0; io<cel.atom[ito].na; ++io)
	{
		double dis=10000.0;

		if(INPUT.system=="hydroxide")
		{
			if(ion==-1) continue;
			if(io!=ion) continue;
				
			for(int i=0; i<36; ++i)
			{				
				double dis = distance(cel.atom[ito].pos[io],cel.atom[itc].pos[i],norm1,norm2,norm3);
				if(dis<INPUT.rcut)
				{
					int index = dis/INPUT.dr;
					assert(index<nmesh);
					trmap[itzone][index]+=1.0;
				}
			}

		//	if(dis2<INPUT.rcut)
		//	{
		//		ofs_running << setw(10) << igeo << setw(10) << io << setw(15) << dis << setw(10) << indexC << setw(15) << dis2 << endl;
		//	}
		}
	}

	delete[] water;

	return;
}
