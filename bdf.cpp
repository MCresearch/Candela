#include "cellFile.h"
#include "input.h"
#include "bdf.h"
#include "math.h"


// general subroutine.
void BDF::Routine()
{
	TITLE("BDF","Routine");

	cout << " Cal bond-angle distribution functions." << endl;

	cal();

	return;
}


// general subroutine.
void BDF::cal()
{
	TITLE("BDF","cal");
	
	//----------------------------------------------------------
	// how many adjacent atoms for each atom we want to mention
	// in this Bond-Angle distribution function.
	//----------------------------------------------------------
	const int bdf_nadj = INPUT.bdf_nadj;
	assert(bdf_nadj > 1 && bdf_nadj < 100);

	//--------------------------------------------
	// information concerning plot of cos(theta).
	//--------------------------------------------

	double dv = INPUT.bdf_dtheta;
	assert( dv > 0.0 );
	int npoints = 180.0/dv + 1;
	cout << " delta value of cos is = " << dv;
	cout << " number of points between -1 and 1 : " << npoints << endl;
	int *bond_df = new int[npoints]; // bond angle distribution function
	for(int i=0; i<npoints; ++i) bond_df[i] = 0;


	// mohan add 2014-09-19
	int nbatype = INPUT.ntype * INPUT.ntype * INPUT.ntype;
	cout << " number of bond angle types :" << nbatype << endl;
	int **bond_detail = new int*[npoints];
	for(int i=0; i<npoints; ++i)
	{
		bond_detail[i] = new int[nbatype];
		for(int j=0; j<nbatype; ++j)
		{
			bond_detail[i][j] = 0;
		}
	}


	// mohan add 2014-10-22
	int* bond_type = new int[INPUT.natom];

	ofstream ofsmovie;
	if(RANK==0 and INPUT.bdf_movie>0)
	{
		ofsmovie.open("BDF_Movie.xyz");
	}

	//--------------------------------------------------
	// calculate the Bond-Angle distribution functions. 
	//--------------------------------------------------
	assert(INPUT.geo_interval>0);
    int count_geometry_number=0;
	
	cout << " geo_1 is " << INPUT.geo_1 << endl;
	cout << " geo_2 is " << INPUT.geo_2 << endl;

    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; igeo=igeo+INPUT.geo_interval)
    {
//		if(igeo%INPUT.geo_interval!=0) continue;

        CellFile cel;

		// get the file name
        stringstream ss; ss << igeo;
        cel.file_name = ss.str();

        // cel : read in geometry file
        if( !CellFile::ReadGeometry( cel ) ) continue;
		if(igeo<INPUT.geo_ignore)//qianrui add
        {
                        cout<<"ignore:"<<igeo<<endl;
                        cel.clean();
                        continue;
        }
        ++count_geometry_number;

		//------------------------------------------
		// mohan add 2014-10-22
		if (INPUT.bdf_movie>0)
		{
			for(int iat=0; iat<INPUT.natom; ++iat)
			{
				bond_type[iat] = -1;
			}
		}
		//------------------------------------------


        // calculate the Bond-Angle distribution function
		// and save it into bond_df.
		this->get_adjacent_atom_positions(
		    bdf_nadj,	
			cel,
			dv,
			npoints,
			bond_df,
			bond_detail, //mohan add bond_detail 2014-09-19	
			bond_type); // mohan add 2014-10-22

		//------------------------------------------
		// mohan add 2014-10-22
		// still may have buts for parallel version/
		// mohan note 2016-12-26
		if (INPUT.bdf_movie>0 and RANK==0)
		{
			ofsmovie << INPUT.natom << endl;
			ofsmovie << "title" << endl;
			int iat=0;
			for(int it=0; it<INPUT.ntype; ++it)
			{
				for(int ia=0; ia<cel.atom[it].na; ++ia)
				{
					ofsmovie << "bond" << bond_type[iat] << " " << cel.atom[it].pos[ia].x << " " 
						<< cel.atom[it].pos[ia].y << " "
						<< cel.atom[it].pos[ia].z << endl;
					++iat;
				}
			}
		}
	}



	int count_total_bond_angles = 0;
	for(int i=0; i<npoints; ++i)
	{
		count_total_bond_angles += bond_df[i];
	}
#ifdef __MPI
	// data 1
	int tmp = count_total_bond_angles;
	count_total_bond_angles=0;
	MPI_Allreduce(&tmp, &count_total_bond_angles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// data 2
	int* tmp2 = new int[npoints];
	for(int i=0;i<npoints;++i){ tmp2[i]=bond_df[i]; bond_df[i]=0;}
    MPI_Allreduce(tmp2, bond_df, npoints, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	delete[] tmp2;
	// data 3
	int **tmp3 = new int*[npoints];
	for(int i=0; i<npoints; ++i)
	{
		tmp3[i] = new int[nbatype];
		for(int j=0; j<nbatype; ++j)
		{
			tmp3[i][j] = bond_detail[i][j];
		}
	}
	for(int i=0; i<npoints; ++i)
	{
		MPI_Allreduce(tmp3[i], bond_detail[i], nbatype, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}
	for(int i=0; i<npoints; ++i)
	{
		delete[] tmp3[i];
	}
	delete[] tmp3;
#endif


	// divided by 'count_total_bond+angles' to 
	// get percent of total bond angles.
	
	// print the Bond-Angle distribution functions.
	if(RANK==0)
	{
		ofstream ofs("BDF.txt");
		ofs << "angle total" ;
		for(int i=0; i<INPUT.ntype; ++i)
		{
			for(int j=0; j<INPUT.ntype; ++j)
			{
				for(int k=0; k<INPUT.ntype; ++k)
				{
					stringstream ijk;
					ijk << i << j << k;
					ofs << " " << ijk.str();
				}
			}
		}
		ofs << endl;
		for(int i=0; i<npoints; ++i)
		{
			ofs << i*dv 
				//		<< " " << bond_df[i]
				<< " " << bond_df[i]/(double)count_total_bond_angles;

			for(int j=0; j<nbatype; ++j)
			{
				ofs << " " << bond_detail[i][j]/(double)count_total_bond_angles;
			}

			ofs << " " << endl;
		}
		ofs.close();	
	}
	
	if(count_geometry_number>0)
	{
		cout << " Final BondAngles/Geometry is " << count_total_bond_angles/count_geometry_number << endl; 
	}

	// close movie
	if(RANK==0 and INPUT.bdf_movie>0)
	{
		ofsmovie.close();
	}

	delete[] bond_df;
	delete[] bond_type;
	for(int i=0; i<npoints; ++i)
	{
		delete[] bond_detail[i];
	}
	delete[] bond_detail;

	return;
}


void BDF::get_adjacent_atom_positions(
	const int &bdf_nadj,
	const Cell &cel,
	const double &dv,
	const int npoints,
	int* bond_df,
	int** bond_detail, // mohan add 2014-09-19
	int* bond_type) // mohan add 2014-10-22
{
	TITLE("BDF","get_adjacent_atom_positions");

	//----------------------------------------
	// information concerning adjacent atoms.
	//----------------------------------------
	//cout << " adjacent atoms = " << bdf_nadj << endl;

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
	const int ncell_1 = 1;
	const int ncell_2 = 1;
	const int ncell_3 = 1;

	// (3) calculate the distance between atoms.
	double x2, y2, z2; // atom positions for atom 2.
	double dx, dy, dz; // delta x,y,z between atom1, atom2.

    int iat=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			
#ifdef __MPI
			if( iat%NPROC!=RANK ) {iat++; continue;}	
			ofs_running << setw(10) << "Atom" << setw(10) << iat << endl;
#endif

			if(INPUT.ion_analysis)
			{
				if(!is_ion(it,ia)) continue;
			}

//			cout << cel.atom[it].pos[ia].x << " "
//			<< cel.atom[it].pos[ia].y << " "
  //  		<< cel.atom[it].pos[ia].z << endl;
//  			cout << " For atom " << it << " " << ia << endl;

			// initialize information for atom 1. 
			for(int iadj=0; iadj<bdf_nadj; ++iadj)
			{
				posx[iadj] = 1.0e8;
				posy[iadj] = 1.0e8;
				posz[iadj] = 1.0e8;
				dist2[iadj] = 1.0e8;
				itindex[iadj] = 1e8;
				iaindex[iadj] = 1e8;
			}

			// !! IMPORTANT !!
			// start searching from first atom!
			// not from (it,ia) because
			// the Bond-Angle distribution functions
			// depends on three atoms, so..
            for(int it2=0; it2<INPUT.ntype; ++it2)
            {
                for(int ia2=0; ia2<cel.atom[it2].na; ++ia2)
                {
					// each atom can be interacted with itself.
					// mohan add 2014-12-11
            //      if(it==it2 && ia==ia2) continue;

                    // search for 27 cells.
					double dis=0.0;
					double dis_saved = 1.0e6;
					double x2_saved, y2_saved, z2_saved;
					double it_saved;
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
                                dis = dx*dx + dy*dy + dz*dz;

								// if this atom2 is an adjacent atom,
								// put its coordinate into posx,posy,posz,
								// also put its distance to atom1 into distance.
								for(int iadj=0; iadj<bdf_nadj; ++iadj)
								{
									if( dis > dist2[iadj])
									{
										continue;
									}
									else if( dis <= dist2[iadj] && dis > 0.0)
									{
										shift2back(bdf_nadj,iadj,dist2,posx,posy,posz,itindex,iaindex);
										// we need to record the 'distance' because
										// we would like to further used this array
										// to include more 'close' atoms.
										dist2[iadj] = dis;
										posx[iadj] = x2;
										posy[iadj] = y2;
										posz[iadj] = z2;
										itindex[iadj] = it2;
										iaindex[iadj] = ia2;
										break;
									}
								}// end iadj

							}// end k
						}// end j
					}// end i



				}// end ia2
			}// end it2

			// FOR TEST
			/*
			cout << " after storage " << endl;
			for(int iadj=0; iadj<bdf_nadj; ++iadj)
			{
				cout << setw(5) << iadj <<
					setw(12) << posx[iadj] << 
					setw(12) << posy[iadj] <<
					setw(12) << posz[iadj] <<
					setw(12) << dist2[iadj] << 
					setw(5) << itindex[iadj] << 
					setw(5) << iaindex[iadj] << endl;
			}
			*/


            // FOR TESTS
/*            
			cout << " Atom(type,index) " << it+1 << " " << ia+1 << endl; 
			cout << setw(15) << cel.atom[it].pos[ia].x << setw(15) << cel.atom[it].pos[ia].y 
			<< setw(15) << cel.atom[it].pos[ia].z << endl;
 			for(int iadj=0; iadj<bdf_nadj; ++iadj)
			{
				cout << setw(15) << posx[iadj] << setw(15) << posy[iadj] << setw(15) << posz[iadj] 
				<< setw(15) << sqrt(dist2[iadj]) 
				<< setw(15) << itindex[iadj]
				<< endl;
			}
			*/
			

			this->bond_angle(
			    INPUT.ntype,
			    it,
				iat,
				itindex,
				cel.atom[it].pos[ia].x, // cartesian coordinates of atom1
				cel.atom[it].pos[ia].y,
				cel.atom[it].pos[ia].z,
				bdf_nadj,
				dv,
				posx,posy,posz, // cartesian coordinates of atom2
				bond_df,
				bond_detail,
				bond_type);

            ++iat;
        }// end ia1
	}//end it1

	// --> FINISH <--
	delete[] posx;
	delete[] posy;
	delete[] posz;
	delete[] dist2;
	delete[] itindex;
	delete[] iaindex;

	return;
}


void BDF::bond_angle(
	const int &ntype,
	const int &it1,
	const int &iat,
	int *itindex,
	const double &x1,
	const double &y1,
	const double &z1,
	const int &bdf_nadj,
	const double &dv,
	const double* posx,
	const double* posy,
	const double* posz,
	int* bond_df,
	int** bond_detail,
	int* bond_type)
{
	double a,b,c;
	// adj1 and adj2 are two different adjacent atoms
	// of atom (x1,y1,z1).
	for(int adj1=0; adj1<bdf_nadj; ++adj1)
	{
		for(int adj2=adj1+1; adj2<bdf_nadj; ++adj2)
		{
			// norm of atom1-adj1, atom1-adj2, adj2-adj2 atom pairs.
			a = cal_norm(x1,y1,z1,posx[adj1],posy[adj1],posz[adj1]);
			b = cal_norm(x1,y1,z1,posx[adj2],posy[adj2],posz[adj2]);
			c = cal_norm(posx[adj1],posy[adj1],posz[adj1],posx[adj2],posy[adj2],posz[adj2]);
//			cout << "a,b,c=" << a << " " << b << " " << c << endl;
			
			assert(c!=0.0);
			double cos_theta = cal_angle(a,b,c);	
			if(cos_theta>1.0)
			{
				cos_theta=1.0;
			}
			if(cos_theta<-1.0)
			{
				cos_theta=-1.0;
			}
			const double degree = std::acos(cos_theta) * 180.0 / 3.1415926;
			const int index = degree / dv;
			bond_df[index]++;

//			cout << " degree=" << degree << endl;

			int itype = it1*ntype*ntype + itindex[adj1]*ntype + itindex[adj2];
			assert(itype<8); //for test only
			bond_detail[index][itype]++;

			if(INPUT.bdf_movie>0) 
			{
				bond_type[iat]=itype;
			}
		}
	}
	return;
}

double BDF::cal_angle(
	const double &a,
	const double &b,
	const double &c)
{
	return (a*a+b*b-c*c)/(2.0*a*b);
}

double BDF::cal_norm(
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


void BDF::shift2back(const int &bdf_nadj,
	const int &ip_start,
	double* dist2,
	double* posx,
	double* posy,
	double* posz,
	int* itindex,
	int* iaindex)
{
	// if ip_start=0, ip=0,
	//cout << " ip_start=" << ip_start << endl;
	for(int ip=bdf_nadj-2; ip>=ip_start; --ip)
	{
		posx[ip+1] = posx[ip];
		posy[ip+1] = posy[ip];
		posz[ip+1] = posz[ip];
		dist2[ip+1] = dist2[ip];
		itindex[ip+1] = itindex[ip]; // mohan add 2014-09-19
		iaindex[ip+1] = iaindex[ip]; // mohan add 2014-12-14
	}
	return;
}

bool BDF::is_ion(const int &it, const int &ia)
{

}
