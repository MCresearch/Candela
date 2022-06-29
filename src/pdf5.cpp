#include "cellFile.h"
#include "input.h"
#include "pdf5.h"
#include "math.h"
#include "HBs.h"
#include "Honeycutt.h"

void PDF5::Routine()
{
	TITLE("PDF5","Routine");

	cal();

	return;
}


void PDF5::cal()
{
	TITLE("PDF5","cal");

	cout << "Calculate the radial distribution functions for the first few neighbours." << endl;

	// (1) set the basic parameters.
	// delta r in real space.
	this->dr = INPUT.dr;
	assert(dr>0.0);

	// radius cutoff in real space, usually choose a/2,
	// where a is the lattice constant.
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);

	// number of radial mesh.
	this->nmesh = int(rcut / dr) +  1;

	// number of neighbours
	this->nn = INPUT.neighbours;
	assert(nn>0);
	cout << "number of neighbours " << nn << endl;

	// pair distribution functions.
	double** gr = new double*[nn];
	for(int i=0; i<nn; ++i)
	{
		gr[i] = new double[nmesh];
		for(int j=0; j<nmesh; ++j) gr[i][j] = 0.0;
	}
	

	ofs_running << " dr = " << dr << " Angstrom" << endl;
	ofs_running << " rcut = " << rcut << " Angstrom" << endl;
	ofs_running << " nmesh = " << nmesh << endl;

	// ionic density 
	assert(INPUT.geo_interval>0);
	this->count_geometry_number=0;

	cout << INPUT.geo_1 << " " << INPUT.geo_2 << endl;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)  
	{
//		cout << "igeo=" << igeo << endl;
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;

		// calculate the ionic density
		this->rho_ion = INPUT.natom / cel.volume;

		if(count_geometry_number==0 && RANK==0 && cel.volume>0)
		{
			ofs_running << " Volume of the input cell = " << cel.volume << " A^3" << " rho(64H2O)=" << 64*18*1.6605/cel.volume << endl;
			ofs_running << " Average ion density = " << rho_ion << endl;
			// Fermi vector
			double kf = pow(3*PI*PI*rho_ion,1.0/3.0);
	        ofs_running << " Fermi vector = " << kf << endl; 
			ofs_running << " PI=" << PI << endl;

	        assert(kf > 0.0);
		}

		// (2)
		int option = 1;
		this->periodic_pairs( cel, gr, option);
		cout << "igeo=" << igeo << endl;

        cel.clean();
	}



    assert(count_geometry_number>0);
    ofs_running << " count_geometry_number = " << count_geometry_number << endl;

    if(count_geometry_number>0)
    {
        for(int ir=0; ir<nmesh; ++ir)
        {
			for(int i=0; i<nn; ++i)
			{
            	gr[i][ir] /= count_geometry_number;
			}
        }
    }

	// output the final pair distribution function
	if(RANK==0)
	{
		ofstream ofs(INPUT.geo_out.c_str());
		// we output the pair correlation function and static structrue factor
		for(int j=0; j<nmesh-1; ++j)
		{
			ofs << (j+0.5)*dr;
			for(int i=0; i<nn; ++i)
			{
				ofs << " " << gr[i][j];
			}
			ofs << endl;
		}
		ofs.close();
	}


	// --> clean <--
	for(int i=0; i<nn; ++i)
	{
		delete[] gr[i];
	}

	delete[] gr;

	
	return;
}


void PDF5::periodic_pairs( const Cell &cel, double **gr, const int &option)
{	
	// (1) calculate the norm of each lattice vectors.
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();

	// (2) calculate how many more cells need to be convered. 
	assert(rcut>0);
	const int ncell_1 = int(rcut/norm1)+1;
	const int ncell_2 = int(rcut/norm2)+1;
	const int ncell_3 = int(rcut/norm3)+1;
	ofs_running << " ncell is " << ncell_1 << " " << ncell_2 << " " << ncell_3 << endl;

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
		else if(cel.atom[it].id=="H") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
		else if(cel.atom[it].id=="Na") itna=it;
		//cout << cel.atom[it].id << endl;
	}
//	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
//	if(INPUT.ntype==3){ assert(itc>=0); }


	// record the distances to the nearest neighbours
	double* nn_dis = new double[this->nn];
	int* nn_index = new int[this->nn];

	


	Water *water;

	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
//		cout << "setup Water." << endl;
		water = new Water[cel.atom[ito].na];
		Water::nions = 0;
		HBs::setup_water(cel, water);
	}	

	//-----------------------------------------------------------------------
	

	// (3) calculate the distance between atoms.
	double x2, y2, z2; // atom positions for atom 2.
	double dx0, dy0, dz0; // difference of coordinates between atom 1 and atom2
	double dx, dy, dz; // delta x,y,z between atom1, atom2.
	double dis;
	int which;
	bool within_distance;
	int ia_save = -1;
	int iat=-1;
	for(int it=0; it<INPUT.ntype; ++it)
	{

		if(INPUT.ele1!="none" and cel.atom[it].id != INPUT.ele1) continue;

		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			++iat;

			for(int inn=0; inn<nn; ++inn)
			{
				nn_dis[inn] = 1000.0;
				nn_index[inn] = -1;
			}

			int iat2=-1;
			for(int it2=0; it2<INPUT.ntype; ++it2)
			{
				if(INPUT.ele2!="none" and cel.atom[it2].id != INPUT.ele2) continue;

				for(int ia2=0; ia2<cel.atom[it2].na; ++ia2)
				{
					++iat2;

					//---------------------------------------------------
					// criteria for different systems
					if(it==it2 && ia==ia2) continue;

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


					dis = 1000.0;
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
								double tmp_dis = sqrt( dx*dx + dy*dy + dz*dz );
								if(tmp_dis<dis) dis=tmp_dis;
							} // end k
						} // end j
					} // end i

			//		cout << dis << endl;


					for(int inn=0; inn<nn; ++inn)
					{
						if( dis < nn_dis[inn] )
						{
							for(int jnn=nn-1; jnn>inn; --jnn)
							{
								nn_dis[jnn]=nn_dis[jnn-1];	
								nn_index[jnn]=nn_index[jnn-1];
							}
							nn_dis[inn]=dis;
							nn_index[inn]=ia2; // index of oxygens
							break;
						}
					}// inn 


			//		cout << "update" << endl;
			//		for(int inn=0; inn<nn; ++inn)
			//		{
			//			cout << nn_dis[inn] << endl;
			//		}

				}// ia2 
			}// it2


			// mohan added 2018-04-11


			// put the distances into the corresponding arrays
			if(INPUT.func==1)
			{
				for(int inn=0; inn<nn; ++inn)
				{
					//	cout << nn_dis[inn] << endl;
					if(nn_dis[inn] <= rcut)
					{
						int which = int(nn_dis[inn]/dr);
						gr[inn][which] += 1.00;
					}
				}
			}
			// test
			else if( (INPUT.func==2 or INPUT.func==3) and nn==5)
			{
				bool hydrogen_bonded=false;
				int ia5 = nn_index[4];

				// the adjacent four atoms should not form
				// any hydrogen bonds with ia5
				for(int inn=0; inn<4; ++inn)
				{
					int ia3=nn_index[inn];
					for(int ii=0; ii<water[ia3].naccept; ++ii)
					{
						if(water[ia3].acceptO[ii]==ia5)
						{
							hydrogen_bonded=true;
							break;
						}
					}
					for(int ii=0; ii<water[ia3].ndonate; ++ii)
					{
						if(water[ia3].donateO[ii]==ia5)
						{
							hydrogen_bonded=true;
							break;
						}
					}
				}

				// ia should not form hydrogen bonds with ia5
				for(int ii=0; ii<water[ia].naccept; ++ii)
				{
					if(water[ia].acceptO[ii]==ia5)
					{
						hydrogen_bonded=true;
						break;
					}
				}
				for(int ii=0; ii<water[ia].ndonate; ++ii)
				{
					if(water[ia].donateO[ii]==ia5)
					{
						hydrogen_bonded=true;
						break;
					}
				}

				// decide whether or not to change gr
				if(INPUT.func==2)
				{
					if(!hydrogen_bonded)
					{
						for(int inn=0; inn<nn; ++inn)
						{
							//	cout << nn_dis[inn] << endl;
							if(nn_dis[inn] <= rcut)
							{
								int which = int(nn_dis[inn]/dr);
								gr[inn][which] += 1.00;
							}
						}
					}
				}
				else if(INPUT.func==3)
				{
					if(hydrogen_bonded)
					{
						for(int inn=0; inn<nn; ++inn)
						{
							//	cout << nn_dis[inn] << endl;
							if(nn_dis[inn] <= rcut)
							{
								int which = int(nn_dis[inn]/dr);
								gr[inn][which] += 1.00;
							}
						}
					}
				}
			}
			//cout << endl;
		}// ia
	}// it

	
	// delete arrays, classes, clean
	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		delete[] water;
	}

	this->count_geometry_number++;

	// clean up
	delete[] nn_dis;
	delete[] nn_index;

	return;
}

