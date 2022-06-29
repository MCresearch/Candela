#include "cellFile.h"
#include "input.h"
#include "bdf_rcut.h"
#include "math.h"
#include "HBs.h"

// general subroutine
void BDF_Rcut::Routine()
{
	TITLE("BDF_Rcut","Routine");

	cout << " Calculate the bond-angle distribution functions" << endl;

	cal();

	return;
}


// general subroutine
void BDF_Rcut::cal()
{
	TITLE("BDF_Rcut","cal");
	
	//--------------------------------------------
	// Setup of bond angle distribution functions 
	//--------------------------------------------

	// original bond angle distribution function
	double dv = INPUT.bdf_dtheta;
	assert( dv > 0.0 );
	int npoints = 180.0/dv + 1;
	cout << " Interval of angles of Cos is " << dv;
	cout << " Number of points between -1 (180 degrees) and 1 (0 degree): " << npoints << endl;
	double *bond_df = new double[npoints]; // bond angle distribution function
	for(int i=0; i<npoints; ++i) bond_df[i] = 0.0;

	// mohan added 2014-09-19
	// bond angle distribution function for number of types
	int nbatype = INPUT.ntype * INPUT.ntype * INPUT.ntype;
	cout << " number of bond angle types :" << nbatype << endl;
	double **bond_detail = new double*[npoints];
	for(int i=0; i<npoints; ++i)
	{
		bond_detail[i] = new double[nbatype];
		for(int j=0; j<nbatype; ++j)
		{
			bond_detail[i][j] = 0;
		}
	}

	// mohan added 2017-05-15
	ofstream ofs0("r_theta.dat");	
	assert(INPUT.nx>0);
	assert(INPUT.ny>0);
	int nx = INPUT.nx;
	int ny = INPUT.ny;
	this->coord_xy = new double*[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->coord_xy[ix] = new double[ny]();
	}

	// number of radial mesh.
	assert(INPUT.dr>0);
	assert(INPUT.rcut>0);
	this->nmesh = int(INPUT.rcut / INPUT.dr) +  1;
	this->ccc = new double[nmesh]();
	this->aabb = new double*[nmesh];
	for(int i=0; i<nmesh; ++i)
	{
		this->aabb[i] = new double[nmesh]();
	}
		

	// initialization
	this->tot_adj = 0;
	this->tot_count = 0;

	//--------------------------------------------------
	// calculate the Bond-Angle distribution functions. 
	//--------------------------------------------------
	assert(INPUT.geo_interval>0);
    int count_geometry_number=0;
	
	cout << " geo file starts from " << INPUT.geo_1 << endl;
	cout << " geo file ends at " << INPUT.geo_2 << endl;

    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
    {
        CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		// get the file name
        stringstream ss; ss << igeo;
        cel.file_name = ss.str();


        // cel : read in geometry file
        if( !CellFile::ReadGeometry( cel ) ) continue;
        ++count_geometry_number;

	if(cel.read_and_used==false) 
	{
		cel.clean();
		continue;
	}

	cout << " geo=" << igeo << endl; 


        // calculate the Bond-Angle distribution function
		// and save it into bond_df.
	this->get_adjacent_atom_positions(
		cel,
		dv,
		npoints,
		bond_df,
		bond_detail); //mohan add bond_detail 2014-09-19	
	cel.clean();
    }

	// print the Bond-Angle distribution functions.
	ofstream ofs("BDF_rcut.dat");
	int count_total_bond_angles = 0;
	for(int i=0; i<npoints; ++i)
	{
		count_total_bond_angles += bond_df[i];
	}
	cout << "count_total_bond_angles=" << count_total_bond_angles << endl;
	ofs_running << "count_total_bond_angles=" << count_total_bond_angles << endl;

	// divided by 'count_total_bond+angles' to 
	// get percent of total bond angles.
	
	ofs << "angle total" ;
	for(int i=0; i<INPUT.ntype; ++i)
	{
		for(int j=0; j<INPUT.ntype; ++j)
		{
			for(int k=0; k<INPUT.ntype; ++k)
			{
				stringstream ijk;
				ijk << i+1 << j+1 << k+1;
				ofs << " " << ijk.str();
			}
		}
	}
	ofs << endl;

	assert(count_total_bond_angles>0);


	// perform normalization
	double sum=0.0;
	for(int i=0; i<npoints; ++i)
	{
//		double cos1=cos((double)i*dv/180*3.1415926);
//		double cos2=cos((double)(i+1)*dv/180*3.1415926);
//		bond_df[i]/=(cos1-cos2);
		sum += bond_df[i];
	}
	sum *= dv;
	cout << "sum for the whole BDF is " << sum << endl;
	ofs_running << "sum for the whole BDF is " << sum << endl;
	assert(sum>0.0);
	for(int i=0; i<npoints; ++i)
	{
		bond_df[i]/=sum;
	}
	sum=0.0;
	for(int i=0; i<npoints; ++i)
	{
		sum += bond_df[i]*dv;
	}
	cout << "After normalization, new sum for the whole BDF is " << sum << endl;


	for(int j=0; j<nbatype; ++j)
	{
		sum=0.0;
		for(int i=0; i<npoints; ++i)
		{
			double sin_theta = std::sin((double)i*dv/180*3.1415926);
			if(INPUT.func==1) // for water
			{
//				double cos1=cos((double)i*dv/180*3.1415926);
//				double cos2=cos((double)(i+1)*dv/180*3.1415926);
//				bond_detail[i][j]/=(cos1-cos2);
				sum += bond_detail[i][j]*sin_theta;
			}
			else if(INPUT.func==2)
			{
				sum += bond_detail[i][j];
			}
		}
		sum *= dv;
//		cout << "sum " << sum << endl;
		if(sum==0) continue;
		for(int i=0; i<npoints; ++i)
		{
			bond_detail[i][j]/=sum;
		}
		sum=0.0;
		for(int i=0; i<npoints; ++i)
		{
			double sin_theta = std::sin((double)i*dv/180*3.1415926);
			if(INPUT.func==1)
			{
				sum += bond_detail[i][j]*dv*sin_theta;
			}
			else if(INPUT.func==2)
			{
				sum += bond_detail[i][j]*dv;
			}
		}
		cout << "sum for " << nbatype << " is " << sum << endl;
		ofs_running << "sum for " << nbatype << " is " << sum << endl;
	}
	
	//--------------------------------------------------------
	// print out data for bond angle distribution function
	//--------------------------------------------------------
	for(int i=0; i<npoints; ++i)
	{
		ofs << i*dv << " " << bond_df[i];
		
		for(int j=0; j<nbatype; ++j)
		{
			ofs << " " << bond_detail[i][j];
		}
	
		ofs << " " << endl;
	}

	ofs.close();	
	
	if(count_geometry_number>0)
	{
		cout << " Final BondAngles/Geometry is " << count_total_bond_angles/count_geometry_number << endl; 

	}

	if(tot_count>0)
	{
		ofs_running << "average number of adjacent atoms: " << tot_adj/(double)tot_count << endl;
		cout << "average number of adjacent atoms: " << tot_adj/(double)tot_count << endl;
	}

	delete[] bond_df;

	//----------------------------------------
	// plot the distribution for ''r-theta''
	//----------------------------------------
	double sum0=0.0;
	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
//			sum0+=coord_xy[ix][iy]*INPUT.dx*INPUT.dy;
			if(coord_xy[ix][iy]>0)
			{
				coord_xy[ix][iy] = -std::log(coord_xy[ix][iy]*INPUT.factor); 
			}
		}
	}

//	if(sum0>0.0)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			for(int ix=0; ix<nx; ++ix)
			{
				//ofs0 << coord_xy[ix][iy]/sum0 << " "; 
				ofs0 << coord_xy[ix][iy] << " "; 
			}
			ofs0 << endl;
		}
	}
	ofs0.close();

	// plot the distance between adj1 and adj2
	ofstream ofsd("d12.dat");
	sum0=0.0;
	for(int i=0; i<nmesh; ++i)
	{
		sum0 += ccc[i]*INPUT.dr;
	}

	if(tot_count>0 and count_geometry_number>0)
	{
		for(int i=0; i<nmesh; ++i)
		{
			ofsd << INPUT.dr*i << " " << ccc[i]/sum0*tot_count/(double)count_geometry_number << endl;
		}
	}
	ofsd.close();

	// plot the aabb
	ofstream ofsab("aabb.dat");
	for(int i=0; i<nmesh; ++i)
	{
		for(int j=0; j<nmesh; ++j)
		{
			ofsab << aabb[i][j] << " ";
		}
		ofsab << endl;
	}
	ofsab.close();

	// clean up
	for(int i=0; i<nmesh; ++i)
	{
		delete[] aabb[i];
	}
	delete[] aabb;
	
	
	for(int ix=0; ix<nx; ++ix)
	{
		delete[] coord_xy[ix];
	}
	delete[] coord_xy;

	return;
}


void BDF_Rcut::get_adjacent_atom_positions(
	const Cell &cel,
	const double &dv,
	const int npoints,
	double* bond_df,
	double** bond_detail) // mohan add 2014-09-19
{
	TITLE("BDF_Rcut","get_adjacent_atom_positions");

	//-----------------------------------------------------------------------
	// if the system is water, we need to do some hydrogen bond analysis
	// get ito, ith, and itc.
	int ito=-1;
	int ith=-1;
	int itc=-1;
	int itcl=-1;
	int itna=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
		else if(cel.atom[it].id=="Na") itna=it;
	}

	Water *water;

	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
//		cout << "setup Water." << endl;
		water = new Water[cel.atom[ito].na];
		Water::nions = 0;
		HBs::setup_water(cel, water);
	}	
	//cout << "water setup" << endl;
	//----------------------------------------
	// information concerning adjacent atoms.
	//----------------------------------------
	int bdf_nadj=100;

	bool* bonded = new bool[bdf_nadj]; // mohan added 2017-02-23
	for(int i=0; i<bdf_nadj; ++i) bonded[i]=false;
	double* posx = new double[bdf_nadj];
	double* posy = new double[bdf_nadj];
	double* posz = new double[bdf_nadj];
	int* itindex = new int[bdf_nadj];
	int* iaindex = new int[bdf_nadj];
	double* dist2 = new double[bdf_nadj]; //record distance^2
	
	// (1) calculate the norm of each lattice vectors.
 	double a1 = cel.a1.norm();
	double a2 = cel.a2.norm();
	double a3 = cel.a3.norm();
	//cout << " norms of lattice vectors = " << a1 << " " << a2 << " " << a3 << endl;

	// (2) calculate how many more cells need to be convered.
	const int ncell_1 = 2;
	const int ncell_2 = 2;
	const int ncell_3 = 2;

	// (3) calculate the distance between atoms.
	double x2, y2, z2; // atom positions for atom 2.
	double dx, dy, dz; // delta x,y,z between atom1, atom2.

	// (4)
	bool within_distance=false;
	double dx0, dy0, dz0; 
	double rcut = INPUT.bdf_rcut;
	double rcut2 = INPUT.bdf_rcut * INPUT.bdf_rcut;

    int iat=0;
	double avg222=0;
	int count222=0;

	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(INPUT.ele1!="none" and cel.atom[it].id!=INPUT.ele1) {continue;} // mohan added 2016-09-26
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
		
			// mohan added 2017-04-08
			//if(cel.atom[it].pos_ili[ia] > 2.0 or cel.atom[it].pos_ili[ia] < -2.0) continue; 
			//if(cel.atom[it].pos_ili[ia] > 5 or cel.atom[it].pos_ili[ia] < 2) continue; 
			//if(cel.atom[it].pos_ili[ia] > 10.0 or cel.atom[it].pos_ili[ia] < 5.0) continue; 

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// mohan added 2016-10-07
			// make sure this 'it' type atom is the adjacent atoms of 'bdf_center' within cutoff 
			bool close = false;
			for(int it0=0; it0<INPUT.ntype; ++it0)
			{
				//cout << cel.atom[it0].id << endl;
				if(INPUT.bdf_center!="none" and cel.atom[it0].id!=INPUT.bdf_center){continue;}
				for(int ia0=0; ia0<cel.atom[it0].na; ++ia0)
				{
					double dis = distance(cel.atom[it].pos[ia],cel.atom[it0].pos[ia0],a1,a2,a3);
					if(dis<INPUT.rcut1) 
					{
						close=true;
						break;
					}
				}	
			}
			if(INPUT.bdf_center!="none" and !close) continue;


			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			int iadj=0;

			// mohan added 2016-09-26
			if(INPUT.z0!=-10000 or INPUT.z1!=-10000)
			{
				double z = cel.atom[it].pos[ia].z;
				while(z>=INPUT.celldm3) z-=INPUT.celldm3;
				while(z<0) z+=INPUT.celldm3;
				if(z<INPUT.z0) continue;
				if(z>INPUT.z1) continue;
			}



			if(INPUT.func_c==1)
			{
				for(int it2=0; it2<INPUT.ntype; ++it2)
				{
					if(INPUT.ele2!="none" and cel.atom[it2].id!=INPUT.ele2) continue; // mohan added 2016-09-26
					for(int ia2=0; ia2<cel.atom[it2].na; ++ia2)
					{
						if( it==it2 && ia==ia2) continue;

						// check x coordinates
						within_distance = false;
						dx0 = cel.atom[it2].pos[ia2].x - cel.atom[it].pos[ia].x;
						for(int i=-ncell_1; i<=ncell_1; ++i)
						{
							dx = abs( dx0 + i*cel.a1.x );
							if( dx < rcut )
							{
								within_distance = true;
							}
						}
						if(within_distance == false) continue;

						// check y coordinates
						within_distance = false;
						dy0 = cel.atom[it2].pos[ia2].y - cel.atom[it].pos[ia].y;
						for(int i=-ncell_2; i<=ncell_2; ++i)
						{
							dy = abs( dy0 + i*cel.a2.y );
							if( dy < rcut )
							{
								within_distance = true;
							}
						}
						if(within_distance == false) continue;

						// check z coordinates
						within_distance = false;
						dz0 = cel.atom[it2].pos[ia2].z - cel.atom[it].pos[ia].z;
						for(int i=-ncell_3; i<=ncell_3; ++i)
						{
							dz = abs( dz0 + i*cel.a3.z );
							if( dz < rcut )
							{
								within_distance = true;
							}
						}
						if(within_distance == false) continue;

						// search for 27 cells.
						for(int i=-ncell_1; i<=ncell_1; ++i)
						{
							for(int j=-ncell_2; j<=ncell_2; ++j)
							{
								for(int k=-ncell_3; k<=ncell_3; ++k)
								{
									// add cell length
									cel.add_cell_length(it2, ia2, i, j, k, x2, y2, z2);

									// calculate the distance between two atoms |r_1 - r_2|
									dx = cel.atom[it].pos[ia].x - x2;
									dy = cel.atom[it].pos[ia].y - y2;
									dz = cel.atom[it].pos[ia].z - z2;
									double dis = dx*dx + dy*dy + dz*dz;

									if( dis <= rcut2 )
									{
										dist2[iadj] = dis;
										// mohan added new analysis 2017-02-23
										bonded[iadj] = false;
										if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
										{
											for(int idd=0; idd<water[ia].ndonate; ++idd)
											{
												if(ia2==water[ia].donateO[idd])
												{
													bonded[iadj]=true;
												}	
											} 
											for(int idd=0; idd<water[ia].naccept; ++idd)
											{
												if(ia2==water[ia].acceptO[idd])
												{
													bonded[iadj]=true;
												}	
											} 
										}
										posx[iadj] = x2;
										posy[iadj] = y2;
										posz[iadj] = z2;
										itindex[iadj] = it2;
										iaindex[iadj] = ia2;
										++iadj;
										//	cout << " iadj=" << iadj << endl;
									}
								}
							}
						}
					}// end ia2
				}// end it2
			}
			// about bonding pairs of wannier functions
			else if(INPUT.func_c==2)
			{
				// mohan added 2017-04-13
				for(int ib=0; ib<INPUT.nbands; ++ib)
				{
					// distance between Wannier centers and oxygen atoms
					double dx=shortest(cel.atom[it].pos[ia].x, cel.wan_centers[ib].x, INPUT.celldm1);
					double dy=shortest(cel.atom[it].pos[ia].y, cel.wan_centers[ib].y, INPUT.celldm2);
					double dz=shortest(cel.atom[it].pos[ia].z, cel.wan_centers[ib].z, INPUT.celldm3);
					double dis = sqrt(dx*dx+dy*dy+dz*dz);
					//if(dis<0.425)
					if(dis>0.425 and dis<0.8)
					{
						posx[iadj] = cel.atom[it].pos[ia].x-dx;
						posy[iadj] = cel.atom[it].pos[ia].y-dy;
						posz[iadj] = cel.atom[it].pos[ia].z-dz;
						++iadj;
					}
				}
				assert(iadj==2);
			}

//			ofs_running << "Number of adjacent atoms is " << iadj << endl;
			tot_adj += iadj;
			tot_count += 1;	
			avg222 += iadj;
			count222 += 1;

			this->bond_angle(
			    INPUT.ntype,
			    it,
				iat,
				itindex,
				cel.atom[it].pos[ia].x, // cartesian coordinates of atom1
				cel.atom[it].pos[ia].y,
				cel.atom[it].pos[ia].z,
				iadj,
				dv,
				bonded,
				posx,posy,posz, // cartesian coordinates of atom2
				bond_df,
				bond_detail);

            ++iat;
        }// end ia1
	}//end it1
	
	if(count222>0)
	{
		ofs_running << "vol_vs_hb " << cel.volume << " " << avg222/(double)count222 << endl;
	}

	// --> FINISH <--
	delete[] bonded;
	delete[] posx;
	delete[] posy;
	delete[] posz;
	delete[] dist2;
	delete[] itindex;
	delete[] iaindex;

	// delete arrays, classes, clean
	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		delete[] water;
	}

	return;
}


void BDF_Rcut::bond_angle(
	const int &ntype,
	const int &it1,
	const int &iat,
	int *itindex,
	const double &x1,
	const double &y1,
	const double &z1,
	const int &bdf_nadj,
	const double &dv,
	const bool* bonded,
	const double* posx,
	const double* posy,
	const double* posz,
	double* bond_df,
	double** bond_detail)
{
	double a,b,c;
	// adj1 and adj2 are two different adjacent atoms
	// of atom (x1,y1,z1).
	for(int adj1=0; adj1<bdf_nadj; ++adj1)
	{
		for(int adj2=adj1+1; adj2<bdf_nadj; ++adj2)
		{
			bool bonded_go_on=true;
			// default is 1, do nothing
			if(INPUT.func_b==2)
			{
				// both bonded
				if(!bonded[adj1] or !bonded[adj2]) continue;
			}
			if(INPUT.func_b==3)
			{
				// both not bonded 
				if(bonded[adj1] or bonded[adj2]) continue;
			}
			if(INPUT.func_b==4)
			{
				// one bonded, one not bonded 
				if(bonded[adj1] and bonded[adj2]) continue;
				if(!bonded[adj1] and !bonded[adj2]) continue;
			}

			// norm of atom1-adj1, atom1-adj2, adj2-adj2 atom pairs.
			//a = cal_norm(x1,y1,z1,posx[adj1],posy[adj1],posz[adj1]);
			//b = cal_norm(x1,y1,z1,posx[adj2],posy[adj2],posz[adj2]);
			//c = cal_norm(posx[adj1],posy[adj1],posz[adj1],posx[adj2],posy[adj2],posz[adj2]);
			double ax = shortest(x1,posx[adj1],INPUT.celldm1);
			double ay = shortest(y1,posy[adj1],INPUT.celldm2);
			double az = shortest(z1,posz[adj1],INPUT.celldm3);
			a = sqrt(ax*ax + ay*ay + az*az);

			double bx = shortest(x1,posx[adj2],INPUT.celldm1);
			double by = shortest(y1,posy[adj2],INPUT.celldm2);
			double bz = shortest(z1,posz[adj2],INPUT.celldm3);
			b = sqrt(bx*bx + by*by + bz*bz);

			double cx = shortest(posx[adj1],posx[adj2],INPUT.celldm1);
			double cy = shortest(posy[adj1],posy[adj2],INPUT.celldm2);
			double cz = shortest(posz[adj1],posz[adj2],INPUT.celldm3);
			c = sqrt(cx*cx + cy*cy + cz*cz);




//			cout << "a,b,c=" << a << " " << b << " " << c << endl;
			
			assert(c!=0.0);
			//double cos_theta = cal_angle(a,b,c);	
			double cos_theta = (ax*bx+ay*by+az*bz)/a/b;
			if(cos_theta>1.0)
			{
				cos_theta=1.0;
				cout << "warning! cos_theta>1.0 " << cos_theta << endl;
			}
			if(cos_theta<-1.0)
			{
				cos_theta=-1.0;
				cout << "warning! cos_theta<-1.0 " << cos_theta << endl;
			}
			double degree = std::acos(cos_theta) * 180.0 / 3.1415926;
			if(degree>180 and degree<360) degree = 360-degree; // mohan added 2017-03-08


			//mohan added 2017-05-15
			if(INPUT.func_b==3)
			{
				if(degree < 50)
				{
					assert(INPUT.dr>0);
					int index00 = c/INPUT.dr;
					if(index00 < this->nmesh)
					{
						this->ccc[index00]+=1; 
					}

					int index11 = a/INPUT.dr;
					int index22 = b/INPUT.dr;
					if(index11 < nmesh and index22 < nmesh)
					{
						this->aabb[index11][index22]+=1; 
					}
				}
			}


			// a typical snapshot
			/* 
			if(INPUT.func_b==3)
			{
				if(degree>54 and degree<55)
				{
					ofs_running << "center atom index=" << iat+1 << " " << adj1 << " " << adj2 << " " << degree << endl;;
				}
			}
			*/

			const int index = degree / dv;

			bond_df[index]+=1.0;


			// calcualte the r-theta 
			assert(INPUT.dy>0);
			int indexY = (degree-INPUT.y0)/INPUT.dy;
			int indexX = (a-INPUT.x0)/INPUT.dx;
			//cout << indexX << " " << indexY << endl;
			if(indexX<INPUT.nx and indexY<INPUT.ny)
			{
				coord_xy[indexX][indexY]++;	
			}
			indexX = (b-INPUT.x0)/INPUT.dx;
			if(indexX<INPUT.nx and indexY<INPUT.ny)
			{
				coord_xy[indexX][indexY]++;	
			}


			if(INPUT.func_c==1)
			{
				int itype = it1*ntype*ntype + itindex[adj1]*ntype + itindex[adj2];
				assert(itype<27); //for test only
				bond_detail[index][itype]+=1.0;
			}
		}
	}
	return;
}

double BDF_Rcut::cal_angle(
	const double &a,
	const double &b,
	const double &c)
{
	return (a*a+b*b-c*c)/(2.0*a*b);
}

double BDF_Rcut::cal_norm(
	const double &x1,
	const double &y1,
	const double &z1,
	const double &x2,
	const double &y2,
	const double &z2)
{
	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;
	return sqrt( dx*dx + dy*dy + dz*dz );
}

