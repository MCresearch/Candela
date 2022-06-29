#include "cellFile.h"
#include "input.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"
#include "waterwire2.h"

int Wire::n=0;

Wire::Wire()
{
	o1=o2=o3=0;
	h1=h2=0;
	f_active=false;
}

Wire::~Wire()
{

}

Waterwire2::Waterwire2()
{
	avg_d=0.0;
	count_d=0;	
    ss_index=0;
}

Waterwire2::~Waterwire2(){}

void Waterwire2::Routine()
{
	TITLE("Waterwire2","Routine");
	
	cout << "Compute the Properties of Waterwire2" << endl;

	ifs_wire.open("wire.dat");
	if(!ifs_wire)
	{
		cout << "Cannot find the wire.dat" << endl;
	}

	ww = new Wire[50000];

	// OPEN FILES 
	PT.setup_PT();

	assert(INPUT.nx>0);
	assert(INPUT.ny>0);
	int nx = INPUT.nx;
	int ny = INPUT.ny;

	this->coord_xy = new double*[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->coord_xy[ix] = new double[ny]();
	}
	// mohan added 2017-06-10
	this->freeE = new double*[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->freeE[ix] = new double[ny]();
		for(int iy=0; iy<ny; ++iy)
		{
			this->freeE[ix][iy]=1.0e6;
		}
	}
	ifstream ifsf("freeE.dat");
	if(ifsf)
	{
		cout << "read in free energy" << endl;
		cout << "ref_rho=" << INPUT.ref_rho << endl;
		cout << "dz=" << INPUT.dz << endl;
		cout << "reference value " << INPUT.ref_rho + INPUT.dz << endl;
		for(int iy=0; iy<ny; ++iy)
		{
			for(int ix=0; ix<nx; ++ix)
			{
				ifsf >> freeE[ix][iy];
				if(freeE[ix][iy] < INPUT.ref_rho + INPUT.dz)
				{
					cout << setw(5) << ix << setw(5) << iy << setw(10) << freeE[ix][iy] << endl;
				}
			}
		}
	}


	int nxx = INPUT.rcut/INPUT.dr;
	dis_o123 = new double[nxx]();
	avg123 = 0.0;
	count123 = 0;


	// begin reading wire information

	this->read_wire();

	// BEGIN CALCULATING DATA
	this->count_geometry_number=0;
	int ipt=0;
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

		if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;
		search_compression(cel, igeo);
		cel.clean();
	}

	ifs_wire.close();

	//
	ofstream ofs_o123("dis_o123.dat");
	double sum=0.0;
	for(int i=0; i<nxx; ++i)
	{
		sum += dis_o123[i]*INPUT.dr;
	}
	if(sum>0.0)
	{
		for(int i=0; i<nxx; ++i)
		{
			dis_o123[i]/=sum;
		}
	}
	for(int i=0; i<nxx; ++i)
	{
		double free_energy = -std::log(dis_o123[i]*sum*INPUT.factor);
		ofs_o123 << i*INPUT.dr << " " << dis_o123[i] << " " << free_energy << endl;
	}
	ofs_o123.close();
	delete[] dis_o123;
	if(count123>0)
	{
		ofs_running << "avg123 (Angstroms) is " << avg123/(double)count123 << endl; 
	}


	// print out data
	ofstream ofs_result("waterwire.dat");
	ofstream ofs_prob("waterwire_coord_prob.txt"); // renxi 20211212

	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			ofs_prob << coord_xy[ix][iy] << " ";
		}
		ofs_prob << endl;
	}
	ofs_prob.close();
	/*
	double sum = 0;
	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			sum += coord_xy[ix][iy];
		}
	}
	sum *= INPUT.dx * INPUT.dy;
	*/
	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			if(coord_xy[ix][iy]>0)
			{
				coord_xy[ix][iy]=-std::log(coord_xy[ix][iy])*INPUT.factor; // INPUT.factor should be kBT. renxi 20211212
			}
		}
	}

	double minleft=1.0e5;
	double minright=1.0e5;
	double ave0=0.0;
	double ave1=0.0;
	int count0=0;
	int count1=0;
	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			ofs_result << coord_xy[ix][iy] << " "; 

			if(ix<nx/2 and coord_xy[ix][iy]<minleft)
			{
				minleft=coord_xy[ix][iy];
			}
			else if(ix>=nx/2 and coord_xy[ix][iy]<minright)
			{
				minright=coord_xy[ix][iy];
			}

			double xx = INPUT.x0 + INPUT.dx*ix;
			double yy = INPUT.y0 + INPUT.dy*iy;

			if( xx > -0.85 and xx < -0.65 and yy > 5.4 and yy < 5.6)
			{
				ave1 += coord_xy[ix][iy];
				++count1;
			}

			if( xx > -0.1 and xx < 0.1 and yy > 4.9 and yy < 5.1)
			{
				ave0 += coord_xy[ix][iy];
				++count0;
			}
		
		}
		ofs_result << endl;
	}
	ofs_result.close();
	cout << "minimal free energy on left is " << minleft << endl;
	cout << "minimal free energy on right is " << minright << endl;
	cout << "middle part x[-0.1,0.1] y[4.9,5.1] is " << ave0/(double)count0 << endl;
	cout << "left part x[-0.85,-0.65] y[5.4,5.6] is " << ave1/(double)count1 << endl;

	if(count_d>0)
	{
		cout << "avg_d is " << avg_d/(double)count_d << endl;
	}

	for(int iw=0; iw<Wire::n; ++iw)
	{
		ofs_running << setw(10) << iw << " o123 " 
		<< ww[iw].o1+1 << " " 
		<< ww[iw].o2+1 << " "
		<< ww[iw].o3+1 << " " 
		<< ww[iw].h1+1 << " " 
		<< ww[iw].h2+1 << " " 
		<< ww[iw].f_active << " "
		<< endl;
	}


	// clean up
	for(int ix=0; ix<nx; ++ix)
	{
		delete[] coord_xy[ix];
	}
	delete[] coord_xy;
	// free energy
	for(int ix=0; ix<nx; ++ix)
	{
		delete[] freeE[ix];
	}
	delete[] freeE;

	delete[] ww;

	return;
}


void Waterwire2::search_compression(const Cell &cel, const int &igeo)
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

	if(INPUT.func==2)
	{
		//all_ions_involved_water(cel, water, ito, ith);
		all_water(cel, water, ito, ith);
	}
	
	double_jump(cel, water, ito, ith);

	delete[] water;
}


void Waterwire2::read_wire()
{
	while(ifs_wire.good())
	{	
		ifs_wire >> ss_index >> ss_time;
		//ofs_running << ss_index << endl;	

		int nw=0;
		ifs_wire >> nw;
		//	cout << "number of water wires is " << nw << endl; 

		string tmp;
		for(int ir=0; ir<nw; ++ir)
		{
			int io1, io2, io3, ih1, ih2;
			ifs_wire >> tmp >> io1 >> io2 >> io3 >> tmp >> ih1 >> ih2;
			io1-=1; io2-=1; io3-=1; ih1-=1; ih2-=1;
			bool f_exist=false;
			for(int iw=0; iw<Wire::n; ++iw)
			{
				if(io1==ww[iw].o1 and io2==ww[iw].o2 and io3==ww[iw].o3)
				{
					if(ih1==ww[iw].h1 and ih2==ww[iw].h2)
					{ 
						f_exist=true;
					}
				}
			}
			if(!f_exist)
			{
				ww[Wire::n].o1=io1;
				ww[Wire::n].o2=io2;
				ww[Wire::n].o3=io3;
				ww[Wire::n].h1=ih1;
				ww[Wire::n].h2=ih2;
				ww[Wire::n].f_active=true;
				++Wire::n;
				cout << "new waterwire " << Wire::n << endl;
			}
		}

	}

	// eliminate double counting
	for(int ip=1; ip<PT.npt-2; ++ip)
	{
		if(PT.type_pt[ip-1]=="double")
		{
			// indexes for atoms on the waterwire
			int o1 = PT.ion_index[ip-1]-1;
			int o2 = PT.ion_index[ip]-1;
			int o3 = PT.ion_index[ip+1]-1;
			int h1 = PT.indexH[ip]-1;
			int h2 = PT.indexH[ip+1]-1;
			for(int iw=0; iw<Wire::n; ++iw)
			{
				if(ww[iw].o1==o1)
				if(ww[iw].o2==o2)
				if(ww[iw].o3==o3)
				if(ww[iw].h1==h1)
				if(ww[iw].h2==h2)
				ww[iw].f_active=false;
			}
		}
	}


	return;
}

void Waterwire2::all_ions_involved_water(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	// begin computing the distribution
	for(int iw=0; iw<Wire::n; ++iw)
	{
		if(!ww[iw].f_active) continue;
		int o1=ww[iw].o1;
		int o2=ww[iw].o2;
		int o3=ww[iw].o3;
		int h1=ww[iw].h1;
		int h2=ww[iw].h2;

		// we don't need to define HB21 and HB32 here because
		// they only appear after PTs occur
		bool formed_HB12 = false;
		bool formed_HB23 = false;
		if(INPUT.system=="hydronium")
		{
			for(int id=0; id<water[o1].ndonate; ++id) 
					if( water[o1].donateO[id] == o2 and water[o1].donateH[id] == h1) formed_HB12=true;
			for(int id=0; id<water[o2].ndonate; ++id) 
					if( water[o2].donateO[id] == o3 and water[o2].donateH[id] == h2) formed_HB23=true;
		}
		if(INPUT.system=="hydroxide")
		{
			for(int id=0; id<water[o1].naccept; ++id) 
					if( water[o1].acceptO[id] == o2 and water[o1].acceptH[id] == h1) formed_HB12=true;
			for(int id=0; id<water[o2].naccept; ++id) 
					if( water[o2].acceptO[id] == o3 and water[o2].acceptH[id] == h2) formed_HB23=true;
		}

		if(!formed_HB12 or !formed_HB23) continue; 

		// distances between o1 and o2, o2 and o3
		double dis1 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm());
		double dis4 = distance(cel.atom[ito].pos[o2], cel.atom[ito].pos[o3], cel.a1.norm(), cel.a2.norm(), cel.a3.norm());

		// distances between o1 and h1, o2 and h1
		double dis2 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
		double dis3 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

		// distances between o2 and h2, o3 and h2
		double dis5 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
		double dis6 = distance(cel.atom[ito].pos[o3], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

		// proton transfer coordinate 1 and 2
		double v1 = dis2 - dis3;
		double v2 = dis5 - dis6;
		double ptcoord = (v1 + v2)/2.0;

		int indexX=0;
		if(INPUT.system=="hydronium")
		{
			indexX = (ptcoord-INPUT.x0)/INPUT.dx;
		}
		if(INPUT.system=="hydroxide")
		{
			indexX = (-ptcoord-INPUT.x0)/INPUT.dx;
		}

		int indexY = (dis1+dis4-INPUT.y0)/INPUT.dy;

		if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0)
		{
			coord_xy[indexX][indexY]+=1.0;

			double r123 = dis1 + dis4;
			if(r123 < INPUT.rcut)
			{
				int lll = int(r123/INPUT.dr);
				this->dis_o123[lll]+=1.0;
				this->avg123+=r123;
				this->count123+=1;
			}
		}

	}	

	return;
}


void Waterwire2::double_jump(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	for(int ip=1; ip<PT.npt-2; ++ip)
	{
		if(PT.type_pt[ip-1]=="double")
		{
			// PT.snapshot_index_pt[ip-1]: time for h1 (first PT) to leave o1 
			// PT.snapshot_index[ip]: time for h1 to arrive o2
			// PT.snapshot_index_pt[ip]: time for h2 to leave o2
			// PT.snapshot_index[ip+1]: time for h2 to arrive o3
			int s1=PT.snapshot_index_pt[ip-1];
			int s2=PT.snapshot_index[ip+1];
			double t1=PT.snapshot_time_pt[ip-1];
			double t2=PT.snapshot_time[ip+1];

			// indexes for atoms on the waterwire
			int o1 = PT.ion_index[ip-1]-1;
			int o2 = PT.ion_index[ip]-1;
			int o3 = PT.ion_index[ip+1]-1;
			int h1 = PT.indexH[ip]-1;
			int h2 = PT.indexH[ip+1]-1;

			// distance between o1 and o2, and o2 and o3
			double dis1 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			double dis4 = distance(cel.atom[ito].pos[o2], cel.atom[ito].pos[o3], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 


			if(h1==-1) continue;
			if(h2==-1) continue;

			assert(o1>=0);
			assert(o2>=0);
			assert(o3>=0);
			assert(h1>=0);
			assert(h2>=0);

			// distances between o1 and h1, o2 and h1
			double dis2 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			double dis3 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

			// distances between o2 and h2, o3 and h2
			double dis5 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			double dis6 = distance(cel.atom[ito].pos[o3], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

			// distance between o1 and o3
			double dis7 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o3], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

			// waterwire length
			double r123 = dis1+dis4;

			// proton transfer coordinate 1 and 2
			double v1 = dis2 - dis3;
			double v2 = dis5 - dis6;
			double ptcoord = (v1 + v2)/2.0;

			bool should_count=true;
			if(INPUT.system=="hydronium")
			{
				// before first jump
				if(cel.snapshot_index < s1 and (v1>0 or v2>0))
				{
					should_count=false;
				}
				// after second jump
				if(cel.snapshot_index > s2 and (v1<0 or v2<0))
				{
					should_count=false;
				}
				// be within 0.5 ps, mohan added 20170429
				if(v1>0 and v2<0)
				{
					if(cel.snapshot_time < t1 or cel.snapshot_time > t2)
					{
						should_count=false;
					}
				}

				bool formed_HB12 = false;
				bool formed_HB21 = false;
				bool formed_HB23 = false;
				bool formed_HB32 = false;
				for(int id=0; id<water[o1].ndonate; ++id) 
					if( water[o1].donateO[id] == o2 and water[o1].donateH[id] == h1) formed_HB12=true;
				for(int id=0; id<water[o2].ndonate; ++id) 
					if( water[o2].donateO[id] == o1 and water[o2].donateH[id] == h1) formed_HB21=true;
				for(int id=0; id<water[o2].ndonate; ++id) 
					if( water[o2].donateO[id] == o3 and water[o2].donateH[id] == h2) formed_HB23=true;
				for(int id=0; id<water[o3].ndonate; ++id) 
					if( water[o3].donateO[id] == o2 and water[o3].donateH[id] == h2) formed_HB32=true;
				if(formed_HB12 and formed_HB23){} 
				else if(formed_HB12 and formed_HB32){}
				else if(formed_HB21 and formed_HB32){}
				else if(formed_HB21 and formed_HB23){}
				else should_count=false;

			}
			if(INPUT.system=="hydroxide")
			{
				if(cel.snapshot_index < s1 and (v1 < 0 or v2 < 0))
				{
					should_count=false;
				}
				if(cel.snapshot_index > s2 and (v1 > 0 or v2 > 0))
				{
					should_count=false;
				}

				// be within 0.5 ps, mohan added 20170429
				if(v1>0 and v2<0)
				{
					if(cel.snapshot_time < t1 or cel.snapshot_time > t2)
					{
						should_count=false;
					}
				}

				// mohan added 2017-03-24
				bool formed_HB12 = false;
				bool formed_HB21 = false;
				bool formed_HB23 = false;
				bool formed_HB32 = false;
				for(int id=0; id<water[o1].naccept; ++id) 
					if( water[o1].acceptO[id] == o2 and water[o1].acceptH[id] == h1) formed_HB12=true;
				for(int id=0; id<water[o2].naccept; ++id) 
					if( water[o2].acceptO[id] == o1 and water[o2].acceptH[id] == h1) formed_HB21=true;
				for(int id=0; id<water[o2].naccept; ++id) 
					if( water[o2].acceptO[id] == o3 and water[o2].acceptH[id] == h2) formed_HB23=true;
				for(int id=0; id<water[o3].naccept; ++id) 
					if( water[o3].acceptO[id] == o2 and water[o3].acceptH[id] == h2) formed_HB32=true;
				if(formed_HB12 and formed_HB23){} 
				else if(formed_HB12 and formed_HB32){}
				else if(formed_HB21 and formed_HB32){}
				else if(formed_HB21 and formed_HB23){}
				else should_count=false;
			}


			int indexX=-1;
			int indexY=-1;

			if(INPUT.system=="hydronium")
			{
				indexX = (ptcoord-INPUT.x0)/INPUT.dx;
			}
			if(INPUT.system=="hydroxide")
			{
				indexX = (-ptcoord-INPUT.x0)/INPUT.dx;
			}
			indexY = (r123-INPUT.y0)/INPUT.dy;

			if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0 and should_count==true)
			{

				double r123 = dis1 + dis4;
				if(r123 < INPUT.rcut)
				{
					if(freeE[indexX][indexY] < INPUT.ref_rho + INPUT.dz) // mohan added 2017-06-10
					{
						int lll = int(r123/INPUT.dr);
						this->dis_o123[lll]+=1.0;
						this->avg123+=r123;
						this->count123+=1;
					}
				}

				if(INPUT.func_b==1)
				{
					coord_xy[indexX][indexY]+=1.0;
				}
				else if(INPUT.func_b==31)
				{
					if(water[o1].naccept==3) coord_xy[indexX][indexY]+=1.0;
					else if(water[o2].naccept==3) coord_xy[indexX][indexY]+=1.0;
					else if(water[o3].naccept==3) coord_xy[indexX][indexY]+=1.0;
				}
				else if(INPUT.func_b==41)
				{
					if(water[o1].naccept>=4) coord_xy[indexX][indexY]+=1.0;
					else if(water[o2].naccept>=4) coord_xy[indexX][indexY]+=1.0;
					else if(water[o3].naccept>=4) coord_xy[indexX][indexY]+=1.0;
				}
				else if(INPUT.func_b==410)
				{
					if(water[o1].naccept>=4 and water[o1].ndonate==0) coord_xy[indexX][indexY]+=1.0;
					else if(water[o2].naccept>=4 and water[o2].ndonate==0) coord_xy[indexX][indexY]+=1.0;
					else if(water[o3].naccept>=4 and water[o3].ndonate==0) coord_xy[indexX][indexY]+=1.0;
				}
			}
		}
	}
	return;
}


void Waterwire2::all_water(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	for(int o1=0; o1<cel.atom[ito].na; ++o1)
	{

		this->avg_d += water[o1].ndonate;
		this->count_d++;

//		if(INPUT.system=="hydronium" and water[o1].nH!=3) continue; 
//		else if(INPUT.system=="hydroxide" and water[o1].nH!=1) continue;

//		if(water[o1].nH!=2) continue;
	
		int nacc2=0;
		if(INPUT.system=="hydronium") nacc2 = water[o1].ndonate;	
		else if(INPUT.system=="hydroxide") nacc2 = water[o1].naccept;

		for(int iacc2=0; iacc2<nacc2; ++iacc2)
		{
			int o2=-1;
			int h1=-1;
			if(INPUT.system=="hydronium")
			{
				o2 = water[o1].donateO[iacc2];
				h1 = water[o1].donateH[iacc2];
			}
			else if(INPUT.system=="hydroxide")
			{
				o2 = water[o1].acceptO[iacc2];
				h1 = water[o1].acceptH[iacc2];
			}
//			if(water[o2].nH!=2) continue;

			double dis1=0.0;
			if(INPUT.system=="hydronium") dis1 = water[o1].donate_disO[iacc2];
			else if(INPUT.system=="hydroxide") dis1 = water[o1].accept_disO[iacc2];
			if(dis1 >= INPUT.rcut) continue;

			int nacc3=0;
			if(INPUT.system=="hydronium") nacc3=water[o2].ndonate;
			else if(INPUT.system=="hydroxide") nacc3=water[o2].naccept;	

			for(int iacc3=0; iacc3<nacc3; ++iacc3)
			{
				int o3=-1;
				int h2=-1;
				if(INPUT.system=="hydronium")
				{
					o3 = water[o2].donateO[iacc3];
					h2 = water[o2].donateH[iacc3];
				}
				else if(INPUT.system=="hydroxide")
				{
					o3 = water[o2].acceptO[iacc3];
					h2 = water[o2].acceptH[iacc3];
				}

//				if(water[o3].nH!=2) continue;
				if(o3==o1) continue;

				double dis4=0.0;

				if(INPUT.system=="hydronium") dis4 = water[o2].donate_disO[iacc3];
				else if(INPUT.system=="hydroxide") dis4 = water[o2].accept_disO[iacc3];

				double r123 = dis1 + dis4;
				int indexY=-1;
				if(r123 < INPUT.rcut)
				{
					indexY = (r123-INPUT.y0)/INPUT.dy;
				}

				// distances between o1 and h1, o2 and h1
				double dis2 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
				double dis3 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				// distances between o2 and h2, o3 and h2
				double dis5 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
				double dis6 = distance(cel.atom[ito].pos[o3], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				// proton transfer coordinate 1 and 2
				double v1 = dis2 - dis3;
				double v2 = dis5 - dis6;
				double ptcoord = 0.0;
				if(INPUT.system=="hydronium") ptcoord = (v1 + v2)/2.0;
				else if(INPUT.system=="hydroxide") ptcoord = -(v1 + v2)/2.0;

				int indexX = (ptcoord-INPUT.x0)/INPUT.dx;

				if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0)
				{
					coord_xy[indexX][indexY]+=1.0;
					double r123 = dis1 + dis4;
					if(r123 < INPUT.rcut)
					{
//						if( ptcoord > -1.5 and ptcoord < 0.0 and r123 > 5.4 and r123 < 5.6) // rest part, mohan add 2017-06-10
						{
							if(freeE[indexX][indexY] < INPUT.ref_rho + INPUT.dz)
							{
								int lll = int(r123/INPUT.dr);
								this->dis_o123[lll]+=1.0;
								this->avg123+=r123;
								this->count123+=1;
							}
						}
					}
				}
			}
		}
	}

	return;
}

