#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "dist.h"
#include "math.h"
#include "water.h"
#include "pdf.h" // to use compute_delta
#include "hyper.h"

Dist::Dist()
{
}

Dist::~Dist(){}

void Dist::Routine()
{
	cout << "Compute 3D distribution of atoms/MLWFs around ions(hydroxide,Cl)" << endl;
	cout << "The format will be charge density format or wave function format that can be read by VESTA." << endl;
	cout << "Dimensions: " << INPUT.u1 << " " << INPUT.u2 << " " << INPUT.u3 << endl;

	upper_natom=0;
	lower_natom=0;

	assert(INPUT.u1>0);
	assert(INPUT.u2>0);
	assert(INPUT.u3>0);
	assert(INPUT.nx>0);
	assert(INPUT.ny>0);

	// initialize 2d array
	this->drx = INPUT.rcut/INPUT.nx;
	this->dry = INPUT.rcut1/INPUT.ny;

	this->dist2D = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->dist2D[ix] = new double[INPUT.ny]();
	}

	// initialize 3d array
	this->dux = INPUT.rcut1*2/INPUT.u1;
	this->duy = INPUT.rcut1*2/INPUT.u2;
	this->duz = INPUT.rcut/INPUT.u3;

	this->dist3D = new double**[INPUT.u1];
	for(int ix=0; ix<INPUT.u1; ++ix)
	{
		this->dist3D[ix] = new double*[INPUT.u2];
		for(int iy=0; iy<INPUT.u2; ++iy)
		{
			this->dist3D[ix][iy] = new double[INPUT.u3]();
		}
	}

	this->count_geometry_number=0;
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

		compute_dist2D(cel);
		++count_geometry_number;
	}

	double sum=0.0;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			sum += dist2D[ix][iy];
		}
	}
	assert(sum>0.0);

	// calculate natom
	double fac=1.0;
	if(INPUT.system=="chlorine")
	{
		assert(count_geometry_number>0);
		upper_natom /= (double)count_geometry_number;
		lower_natom /= (double)count_geometry_number;
		ofs_running << "average upper natom: " << upper_natom << endl;
		ofs_running << "average lower natom: " << lower_natom << endl;
		fac = upper_natom+lower_natom+1.0; // the 1.0 account for the Cl we chose
		ofs_running << "total coordination number: " << fac << endl;
	}


	// print out the data
	ofstream ofs("dist2D.dat");
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			dist2D[ix][iy] = dist2D[ix][iy]/sum/drx/dry*fac;
			ofs << dist2D[ix][iy] << " ";
		}
		ofs << endl;
	}
	ofs.close();

	 


	// print out 3D information
	compute_dist3D();
	ofstream ofs3("dist3D.dat");
	cout << "dimensions of 3D array are: " << INPUT.u1 << " " << INPUT.u2 << " " << INPUT.u3 << endl;
	// for ixf format
	for(int iz=0; iz<INPUT.u3; ++iz)
	{
		for(int iy=0; iy<INPUT.u2; ++iy)
		{
			for(int ix=0; ix<INPUT.u1; ++ix)
			{
				ofs3 << dist3D[ix][iy][iz] << endl;
			}
		}
	}
	ofs3 << " END_DATAGRID_3D" << endl;
	ofs3 << " END_BLOCK_DATAGRID_3D" << endl;
	ofs3.close();



	// delete data
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		delete[] dist2D[ix];
	}
	delete[] dist2D;

	for(int ix=0; ix<INPUT.u1; ++ix)
	{
		for(int iy=0; iy<INPUT.u2; ++iy)
		{
			delete[] dist3D[ix][iy];
		}
		delete[] dist3D[ix];
	}
	delete[] dist3D;


	return;
}

void Dist::compute_dist2D(const Cell &cel)
{
	// get it index for each element;
	int ito=-1;
	int ith=-1;
	int itcl=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
	}
	if(INPUT.ntype>=2){ assert(ito>=0); assert(ith>=0);}

	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

	HBs::setup_water(cel, water);

	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(INPUT.system=="hydroxide" and water[ia].nH==1)
		{
//			cout << "Find ion!!" << endl;
			const int io2 = ia;
			const int ih3 = water[ia].indexH[0]; 

			// all the atoms around OH- are included
			if(INPUT.func==1)
			{
				for(int it=0; it<INPUT.ntype; ++it)
				{
					if(cel.atom[it].id == INPUT.ele1)
					{
						for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
						{
							if(ito==it and ia==ia2) continue;
							double angle0 = HBs::angle(cel, cel.atom[it].pos[ia2], cel.atom[ito].pos[io2], cel.atom[ith].pos[ih3]);
							double dis = distance(cel.atom[ito].pos[ia], cel.atom[it].pos[ia2], 
									INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							locate(angle0, dis);
						}
					}
				}	
			}
			// compute the 2D distribution function only for
			// accepted and donating HBs.
			else if(INPUT.func==2)
			{
				if( (INPUT.nacc==40 and water[ia].naccept>=4) or
						(INPUT.nacc==30 and water[ia].naccept>=3) or
						(INPUT.nacc==5 and water[ia].naccept==5) or
						(INPUT.nacc==4 and water[ia].naccept==4) or
						(INPUT.nacc==3 and water[ia].naccept==3) ) 
				{
					for(int ia2=0; ia2<water[ia].naccept; ++ia2)
					{
						const int io1 = water[ia].acceptO[ia2];
						double angle0 = HBs::angle(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io2], cel.atom[ith].pos[ih3]);
						double dis = water[ia].accept_disO[ia2];
						locate(angle0, dis);
					}

					for(int ia2=0; ia2<water[ia].ndonate; ++ia2)
					{
						const int io1 = water[ia].donateO[ia2];
						double angle0 = HBs::angle(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io2], cel.atom[ith].pos[ih3]);
						double dis = water[ia].donate_disO[ia2];
						locate(angle0, dis);
					}// end donate
				}
			}
			// planarity parameter is another criterion
			else if(INPUT.func==3)
			{
				cout << water[ia].naccept << endl;
				if(INPUT.nacc==4 and water[ia].naccept==4)	
				{
					Vector3<double> da1, da2, da3, n, d;
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

						if(dis>2.5)
						{
							const int io1 = water[ia].acceptO[i4];
							double angle0 = HBs::angle(cel, cel.atom[ito].pos[io1], cel.atom[ito].pos[io2], cel.atom[ith].pos[ih3]);
							double dis = water[ia].accept_disO[i4];
							locate(angle0, dis);
						}
					}
				}
			}
			// calculate the distributions of wannier functions around OH-
			else if(INPUT.func==4)
			{
				for(int ib=0; ib<cel.nbands; ++ib)
				{
					double dis = distance(cel.atom[ito].pos[ia], cel.wan_centers[ib], 
									INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					if(dis<1.0)
					{
						cout << setw(10) << ib+1 << setw(10) << dis << endl;
						double angle0 = HBs::angle(cel, cel.wan_centers[ib], cel.atom[ito].pos[io2], cel.atom[ith].pos[ih3]);
						locate(angle0, dis);
					}
				}
			}
			else
			{
				cout << "check INPUT.func. " << endl;
				exit(0);
			}
		} // end hydroxide
	} // iao


	// this part is for chlorine
	if(INPUT.system=="chlorine")
	{
		assert(itcl>=0);
		// search for all chlorine atoms
		for(int ia=0; ia<cel.atom[itcl].na; ++ia)
		{
			// search for the neaest neighbor "oxygen" of chlorine
			double shortest_dis=1000;
			int shortest_ia=-1;
			for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
			{
				double dis = distance(cel.atom[itcl].pos[ia], cel.atom[ito].pos[ia2], 
							INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				if(dis<shortest_dis)
				{
					shortest_dis = dis;
					shortest_ia = ia2;
				}
			}	
			cout << "shortest_ia = " << shortest_ia << " shortest_dis = " << shortest_dis << endl;

			if(INPUT.func==1) // search chlorine for center water
			{
				double shortest_disH = 1000;
				int shortest_ih = -1;
				for(int ih=0; ih<water[shortest_ia].nH;++ih)
				{
					int index_H = water[shortest_ia].indexH[ih];
					double dis = distance(cel.atom[itcl].pos[ia], cel.atom[ith].pos[index_H], 
							INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					if(dis<shortest_disH)
					{
						shortest_disH = dis;
						shortest_ih = index_H;
					}
				}
				cout << " shortest_ih = " << shortest_ih << endl;

				double angle0 = HBs::angle(cel, cel.atom[itcl].pos[ia], cel.atom[ito].pos[shortest_ia], cel.atom[ith].pos[shortest_ih]);
				locate(angle0, shortest_dis);
			}
			else if(INPUT.func==2) // search oxygen for center chlorine
			{

				for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
				{
					if(ia2==shortest_ia) continue;
					double dis = distance(cel.atom[itcl].pos[ia], cel.atom[ito].pos[ia2], 
							INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					if(dis<INPUT.rcut_clo)
					{
//						cout << "ia2=" << ia2 << endl;
						double angle0 = HBs::angle(cel, cel.atom[ito].pos[ia2], cel.atom[itcl].pos[ia], cel.atom[ito].pos[shortest_ia]); 
						locate(angle0, dis);
					}
				}
			}
			else if(INPUT.func==3) // search oxygen for center chlorine
			{
				double un=0.0;
				double ln=0.0;
				int nn=0;
				for(int ia3=0; ia3<cel.atom[ito].na; ++ia3)
				{
					double dis3 = distance(cel.atom[itcl].pos[ia], cel.atom[ito].pos[ia3], 
							INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					if(dis3>=INPUT.rcut_clo) continue;
					else nn++;
					for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
					{
						if(ia2==ia3) continue;
						double dis2 = distance(cel.atom[itcl].pos[ia], cel.atom[ito].pos[ia2], 
								INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						if(dis2<INPUT.rcut_clo)
						{
							double angle0 = HBs::angle(cel, cel.atom[ito].pos[ia2], cel.atom[itcl].pos[ia], cel.atom[ito].pos[ia3]); 
							locate(angle0, dis2);
							if(angle0<90) un+=1.0;
							else ln+=1.0;
						}
					}
				}
				assert(nn>1);
				un=un/(double)nn;
				ln=ln/(double)nn;
				cout << "un=" << un << " ln=" << ln << " nn=" << nn << endl;
				upper_natom+=un;
				lower_natom+=ln;
			}
		}
	}

	delete[] water;
}


void Dist::locate(const double &angle0, const double &dis)
{

//	if(dis>=4) return;

	double xx = std::cos(angle0/180*3.1415926)*dis;
	double yy = std::sin(angle0/180*3.1415926)*dis;

	int indx = (double)xx/drx+INPUT.nx/2;
	int indy = (double)yy/dry;

	if(indx<INPUT.nx and indy<INPUT.ny and indx>=0 and indy>=0)
	{
//		cout << indx << " " << indy << endl;
		dist2D[indx][indy] += 1;	
	}
	return;
}


void Dist::compute_dist3D()
{
	cout << "compute dist3D." << endl;

	for(int ix=0; ix<INPUT.u1; ++ix)
	{
		double xx = -INPUT.rcut1 + ix*dux;
//		cout << xx << endl;
		for(int iy=0; iy<INPUT.u2; ++iy)
		{
			double yy = -INPUT.rcut1 + iy*dux;
			for(int iz=0; iz<INPUT.u3; ++iz)
			{
				double disxy = sqrt(xx*xx+yy*yy);
				double disz = iz * this->duz;
				
				const int index_x2d = disz/this->drx;
				const int index_y2d = disxy/this->dry;
			
//				cout << dis2 << " " << index_x2d << endl;
//				cout << dis1 << " " << index_y2d << endl;

				if(index_x2d < INPUT.nx and index_y2d < INPUT.ny)
				{
					dist3D[ix][iy][iz] = dist2D[index_x2d][index_y2d];
				}
			}
		}
	}

	return;
}
