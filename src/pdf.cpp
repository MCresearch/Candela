#include "cellFile.h"
#include "input.h"
#include "pdf.h"
#include "math.h"
#include "HBs.h"
#include "Honeycutt.h"

PDF::PDF(){}

PDF::~PDF(){}

void PDF::Routine()
{
	TITLE("PDF","Routine");

	this->cal();

	return;
}


void PDF::cal()
{
	TITLE("PDF","cal");

	cout << "Calculate the radial distribution functions g(r)." << endl;

	// (1) set 'delta r' in real space from input file.
	this->dr = INPUT.dr;
	assert(dr>0.0);

	// radius cutoff in real space, usually a/2,
	// where a is the lattice constant of a cubic box.
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);

	// number of radial mesh grids.
	this->nmesh = int(rcut / dr) +  1;
	
	// radial distribution function.
	double* gr = new double[nmesh];
	double* gr_tmp = new double[nmesh];
	for(int i=0; i<nmesh; ++i) gr[i] = 0;
	for(int i=0; i<nmesh; ++i) gr_tmp[i] = 0;

	ofs_running << " dr = " << dr << " Angstrom" << endl;
	ofs_running << " rcut = " << rcut << " Angstrom" << endl;
	ofs_running << " nmesh = " << nmesh << endl;

	// prepare for HA index
	if(INPUT.HA_pdf==true)
	{
		cout << "prepare for HA index decomposed PDF" << endl;
		cout << this->HA.max_num_adj << endl;
		HA.nadj = new int[INPUT.natom]();
		HA.adj_index = new int*[INPUT.natom];
		for(int i=0; i<INPUT.natom; ++i)
		{
			HA.adj_index[i] = new int[HA.max_num_adj];
			for(int j=0; j<HA.max_num_adj; ++j)
			{
				HA.adj_index[i][j] = -1;
			}
		}
	}
	cout << 12 << endl;
	//ofs_calc_snapshots.open("snapshots_calculated.txt");
	cout << 2 << endl;
	if (INPUT.nPT > 0)
	{
		//cout << INPUT.nPT << endl;
		ifstream ifs_trans("trans.dat");
		
		//cout << 0 << endl;
		this->iindex = new int[INPUT.nPT];
		this->iindex_p = new int[INPUT.nPT];
		this->snapshot_time_pt = new double[INPUT.nPT];
		this->return_jump = new bool[INPUT.nPT];
		//cout << 1 << endl;
		for (int is=0; is<INPUT.nPT; is++)
		{
			return_jump[is] = false;
		}
		int orig_iidx;
		string useless;
		ifs_trans >> useless >> useless >> useless >> useless >> useless >> useless >> useless >> orig_iidx >> useless >> useless >> useless;
		orig_iidx--;
		cout << 2 << endl;
		for(int is = 0; is<INPUT.nPT; is++)
		{
			int ii, ii_p;
			int ss_pt;
			double ss_time_pt;
			//string pt_type;
			ifs_trans >> useless >> useless >> ii >> useless >> useless >> ss_pt >> ss_time_pt >> ii_p >> useless >> useless >> useless;
			this->iindex[is] = ii-1;
			this->iindex_p[is] = ii_p-1;
			this->snapshot_time_pt[is] = ss_time_pt;
			if (is == 0)
			{
				if (ii-1 == orig_iidx) this->return_jump[is] = true;
			}
			else
			{
				if (ii-1 == this->iindex_p[is-1]) this->return_jump[is] = true;
			}
		}
	}


	// ionic density 
	assert(INPUT.geo_interval>0);
	this->count_geometry_number=0;

	cout << INPUT.geo_1 << " " << INPUT.geo_2 << endl;

	// For calculating std of PDFs.
	int ngeo_count = int((INPUT.geo_2-INPUT.geo_ignore+1)/INPUT.geo_interval);

	if (INPUT.pdf_nstd > 0)
	{
		this->npdf = int(ngeo_count/INPUT.pdf_nstd); // Got to calculate so many pdf.
		this->npdf_count = new int[npdf];
		this->multi_pdf = new double*[npdf];
		for (int ipdf = 0; ipdf<this->npdf; ipdf++)
		{
			this->multi_pdf[ipdf] = new double[int(INPUT.rcut/INPUT.dr) + 1];
			this->npdf_count[ipdf] = 0;
			for (int ir=0; ir<int(INPUT.rcut/INPUT.dr) + 1; ir++)
			{
				this->multi_pdf[ipdf][ir] = 0;
			}
		}
	}


#ifdef __MPI
	int geo_start=(INPUT.geo_ignore<INPUT.geo_1?INPUT.geo_1:INPUT.geo_ignore); //qianrui add
    // total number of used configurations.
    int numConfigs=INPUT.geo_2-geo_start+1;//qianrui change geo_1 to geo_start
	// local number of used configurations.
	int numConfigsLoc = numConfigs/NPROC;
	// extra number of configurations.
	int numConfigsExtra = numConfigs%NPROC;
	// initial configuration
	int numConfigsInit=-100;
    if( RANK < numConfigsExtra )
	{
		numConfigsLoc = numConfigsLoc+1;
		numConfigsInit = RANK * numConfigsLoc + geo_start;//qianrui change geo_1 to geo_start
	}
	else if( RANK >= numConfigs )
	{
		numConfigsLoc = 0;
	}
	else
	{
		numConfigsInit = numConfigsExtra + RANK * numConfigsLoc + geo_start;//qianrui change geo_1 to geo_start
	}

	ofs_running << " numConfigsInit:" << numConfigsInit << endl;
	ofs_running << " numConfigsLoc:" << numConfigsLoc << endl;
	ofs_running << " numConfigsExtra:" << numConfigsExtra << endl;

	for(int igeo=INPUT.geo_1; igeo<numConfigsInit+numConfigsLoc; ++igeo)//qianrui change numConfigsInit to geo_1 
#else
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
#endif
	{
		cout << "igeo=" << igeo << endl;
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;
		cout << "Succeeded" << endl;
		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
#ifdef __MPI //qianrui add
		if(igeo<numConfigsInit) cel.read_and_used=false;
#endif
		if(cel.read_and_used==false) 
		{
			cel.clean();//qianrui add
			continue;
		}
		// calculate the ionic density
		this->rho_ion = INPUT.natom / cel.volume;
		if(count_geometry_number==0 && RANK==0 && cel.volume>0)
		{
			// only works for cell consisting of 64 water molecules
			ofs_running << " Volume of the input cell = " << cel.volume 
			<< " A^3" << " rho(64H2O)=" << 64*18*1.6605/cel.volume << endl;
			ofs_running << " Average ion density = " << rho_ion << endl;
			// Fermi vector
			double kf = pow(3*PI*PI*rho_ion,1.0/3.0);
	        ofs_running << " Fermi vector = " << kf << endl; 
			ofs_running << " pi=" << PI << endl;

	        assert(kf > 0.0);
		}

		// add this 2014-05-30
		for(int i=0; i<nmesh; ++i) gr_tmp[i] = 0;


		// for HA
		if(INPUT.HA_pdf==true)
		{
			Honeycutt::setup_nadj(cel, HA.nadj, HA.adj_index, HA.max_num_adj);
		}
		
		// (2)
		int option = 1;
		//cout << option << endl;
		int index_pdf_in = -1;
		if (INPUT.pdf_nstd > 0)
		{
			for (int ipdf=0; ipdf<this->npdf; ipdf++)
			{
				if (this->count_geometry_number%this->npdf == ipdf and INPUT.func_e == 1)
				{
					index_pdf_in = ipdf;
				}
				if (int(this->count_geometry_number/INPUT.pdf_nstd) == ipdf and INPUT.func_e == 2)
				{
					index_pdf_in = ipdf;
				}
			}	
			//assert(index_pdf_in >= 0);
		}
		this->periodic_pairs( cel, gr_tmp, option, index_pdf_in );
		cout << "igeo=" << igeo << endl;

		double sumt = 0.0;
		for(int i=0; i<nmesh; ++i) {sumt += gr_tmp[i];}
		//	cout << "sumt=" << sumt << endl; 
		ofs_running << "count_geometry_number=" << count_geometry_number << endl;


		// (3) calculate the radial distribution function
		//
		// dV, equals to 4 * pi * r^2 * dr
		const double prec = 4.0/3.0*PI;
		for(int i=0; i<nmesh; ++i)
		{
			// volume
			double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));

			if(INPUT.ntype==1 || INPUT.func_c==2 ) // mohan added 2017-06-23
			{
				//qianrui modified 2020-2-18: In fact -1 won't change the result when N is large enough, but will affect it when N is small.
				//gr_tmp[i] = gr_tmp[i] / INPUT.natom / rho_ion / vv;
				gr_tmp[i] = gr_tmp[i] / (INPUT.natom-1) / rho_ion / vv;
			}
			//else if(INPUT.ntype==2 || (INPUT.ntype==3 and INPUT.system=="water"))
			else if(INPUT.ntype==2 || INPUT.ntype==3) // new on 2017-06-23
			{
				int n1=0;
				int n2=0;
				for(int it=0; it<INPUT.ntype; ++it)
				{
					//cout << "id=" << cel.atom[it].id << endl;
					if(cel.atom[it].id == INPUT.ele1) n1=cel.atom[it].na;
					if(cel.atom[it].id == INPUT.ele2) n2=cel.atom[it].na;
				}
				assert(n1>0);
				assert(n2>0);

				if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
				{
					n1=1;
				}
				if (INPUT.system == "water" and INPUT.nPT > 0 and INPUT.ele1 == "O")
				{
					n1 = 1;
				}

				if(INPUT.func_d==1) // normal
				{
					// rho_ion: ionic density
					// vv: volume
			    	if(INPUT.ele1==INPUT.ele2)//qianrui modified 2020-2-18
						gr_tmp[i] = gr_tmp[i] * INPUT.natom / n1 /(n2-1) / rho_ion / vv;
					else
			    		gr_tmp[i] = gr_tmp[i] * INPUT.natom / n1 /n2 / rho_ion / vv;
				}
				else if(INPUT.func_d==2) // rho*g
				{
					// 2017-12-22
					// the formula in the line below has ionic density information inside
					gr_tmp[i] = gr_tmp[i] * INPUT.natom / n1 /n2 / vv;
				}
			}
		}
		//cout << 1 << endl;
		for(int i=0; i<nmesh; ++i)
		{
			gr[i]+=gr_tmp[i];
		}
		cel.clean();
	}
	//ofs_calc_snapshots.close();
    // do average of 3D ssf.

#ifdef __MPI
    int local=count_geometry_number;
	count_geometry_number=0;
    MPI_Allreduce(&local, &count_geometry_number, 1, MPI_INT , MPI_SUM , MPI_COMM_WORLD);
	double *gr_local = new double[nmesh];
	for(int ir=0; ir<nmesh; ++ir)
	{
		gr_local[ir] = gr[ir];
		gr[ir] = 0.0;
	}
	MPI_Allreduce(gr_local, gr, nmesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
#endif

    assert(count_geometry_number>0);
    ofs_running << " count_geometry_number = " << count_geometry_number << endl;

    if(count_geometry_number>0)
    {
		cout << "count_geometry_number = " << count_geometry_number << endl;
        for(int ir=0; ir<nmesh; ++ir)
        {
            gr[ir] /= count_geometry_number;
        }
		if (INPUT.pdf_nstd > 0)
		{
			this->sort_multi_pdf(this->multi_pdf, this->npdf_count);
		}
    }


	double gr_sum = 0;
	
	// output the final pair distribution function
	if(RANK==0)
	{
		ofstream ofs(INPUT.geo_out.c_str());
		// we output the pair correlation function and static structrue factor
		for(int i=0; i<nmesh-1; ++i)
		{
			ofs << (i+0.5)*dr << " " << gr[i];
			gr_sum += 4*PI*(i+0.5)*(i+0.5)*dr*dr*dr * gr[i] ;
			ofs << " " << gr_sum << endl;
		}
		ofs.close();

		if (INPUT.pdf_nstd > 0)
		{
			this->out_multiple_pdf(this->multi_pdf);
		}
	}


#ifdef __MPI
	delete[] gr_local;
#endif


// calculate the SSF via FFT,
	if(INPUT.struf_dgx>0 and INPUT.struf_ng>0)
	{
		// calculate the static structure factor 
		this->dg = INPUT.struf_dgx;
		this->ng = INPUT.struf_ng;

		double *sf = new double[ng];
		for(int ig=0; ig<ng; ++ig) sf[ig] = 0.0;

		this->static_structure_factor(rho_ion, gr, sf);

		// output the static structure factor.
		ofstream ofss("ssf.txt");
		for(int ig=0; ig<ng; ++ig)
		{
			double k = ig*dg;
			if( k >= 2*PI/(rcut*2.0) )
			{
				ofss << ig*dg << " " << sf[ig] << endl;
			}
		}
		ofss.close();
		delete[] sf;
	}
	
	// --> clean <--

	delete[] gr;
	delete[] gr_tmp;

	// for HA
	if(INPUT.HA_pdf==true)
	{
		for(int i=0; i<INPUT.natom; ++i)
		{
			delete[] HA.adj_index[i]; 
		}
		delete[] HA.adj_index;
		delete[] HA.nadj;
	}

	return;
}


// calculate the static structure factor using direct Fourier Transform.
void PDF::static_structure_factor(const double &rho_ion, const double *gr, double *sf) const
{
	cout << " Begin calculate the structure factor" << endl;
	for(int ig=1; ig<ng; ++ig)
	{
		double k = ig*dg;
		//			cout << " k=" << k << endl;
		double v = 0.0;
		for(int ir=0; ir<nmesh-1; ++ir)
		{
			double r = ((double)ir+0.5)*dr;
			double kr = r * k;
			v += r*r*(gr[ir]-1.0)*sin(kr)/kr*dr;
			//v += r*r*gr[ir]*sin(kr)/kr*dr;
		}
		sf[ig]=1+4*PI*rho_ion*v;
	}
	return;
}


void PDF::periodic_pairs( const Cell &cel, double *gr, const int &option, const int &index_pdf)
{	
	// (1) calculate the norm of each lattice vectors.
	//cout << "start" << endl;
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();
	//cout << "norm1 = " << norm1 << endl;
	// (2) calculate how many more cells need to be convered. 
	assert(rcut>0);
	const int ncell_1 = int(rcut/norm1)+1;
	const int ncell_2 = int(rcut/norm2)+1;
	const int ncell_3 = int(rcut/norm3)+1;
	//ofs_running << " ncell is " << ncell_1 << " " << ncell_2 << " " << ncell_3 << endl;

	bool should_count=false;
	if(INPUT.system!="hydronium" and INPUT.system!="hydroxide")
	{
		should_count=true;
	}

	//-----------------------------------------------------------------------
	// if the system is water, we need to do some hydrogen bond analysis
	// so set we need to set ito, ith, and itc.
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
		//cout << cel.atom[it].id << endl;
	}
//	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
//	if(INPUT.ntype==3){ assert(itc>=0); }


	// setup water molecules if necessary
	Water *water;

	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
//		cout << "setup Water." << endl;
		water = new Water[cel.atom[ito].na];
		Water::nions = 0;
		HBs::setup_water(cel, water);
	}	

	int io_of_ion = -1;
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
	int count_num = 0; // for nPT
	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(INPUT.ele1!="none" and cel.atom[it].id != INPUT.ele1) continue;

		//cout << "it = " << it << endl;

		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			++iat;

			if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
			{
				ia_save = ia;
				// criterion 1: oxygen is in ion
				if(!PDF_ADDED::atom_in_ion(cel,it,ia,io_of_ion)) continue;
				assert(io_of_ion>=0);
				if (INPUT.nPT > 0)
				{
					if (!correct_ion_correct_time(cel, water, ito, io_of_ion)) continue;
					cout << "Calculated_ss" << cel.snapshot_time << " " << io_of_ion << endl;
				}
				if(INPUT.func==1)
				{
					// criterion 2
					if(INPUT.delta!=0.0)
					{
						double delta = PDF_ADDED::compute_delta(cel,water,it,ia);
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
				}
				else if(INPUT.func==2)
				{
					// select the atom that ion donates HB to.
//					if(water[ia].ndonate==1)
//					{
//						ia = water[ia].donateO[0];	
//					}
					cout << "warning! stop." << endl;
					exit(0);
				} 
			}//system
			if(INPUT.system == "water" and INPUT.ele1 == "O")
			{
				if (INPUT.nacc>=0)
				{
					if (water[ia].naccept != INPUT.nacc) continue;
				}
				if (INPUT.ndon>=0)
				{
					if (water[ia].ndonate != INPUT.ndon) continue;
				}
				if (INPUT.nPT > 0)
				{
					if (!correct_ion_correct_time(cel, water, ito, ia)) 
					{
						continue;
					}
					cout << "Calculated_ss" << cel.snapshot_time << " " << ia << endl;
					count_num++;
				}
			}


//			if(cel.atom[it].na > 10000) 
//			{
//				ofs_running << "atom " << ia << endl;
//			}
//			ofs_running << cel.atom[it].pos[ia].x << " " 
//			<< cel.atom[it].pos[ia].y << " "
//			<< cel.atom[it].pos[ia].z << endl;

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

					if(INPUT.system=="water")
					{
						// default is 1
						if(INPUT.func_b==2)
						{
							// only select those intact water molecules
							if(water[ia2].naccept!=2 or water[ia2].ndonate!=2) continue;
						}
						else if(INPUT.func_b==3)
						{
							// only select those NOT intact water molecules
							if(water[ia2].naccept==2 and water[ia2].ndonate==2) continue;
						}
					}
					if (INPUT.func_b == 4)
					{
						assert(INPUT.system=="water" or INPUT.system=="hydroxide" or INPUT.system=="hydronium");
						bool in_donate_water = false;
						for (int iwater=0; iwater<water[ia].ndonate; iwater++)
						{
							int water_index = water[ia].donateO[iwater];
							if (INPUT.ele2 == "O")
							{
								if (water_index == ia2){in_donate_water = true;} 
							}
							else if (INPUT.ele2 == "H")
							{
								for (int iH=0; iH<water[water_index].nH; iH++)
								{
									if (ia2 == water[water_index].indexH[iH])
									{
										in_donate_water = true;
									}
								}
							}
						}
						if(!in_donate_water)
						{
							continue;
						}
					} // renxi 20210607

					else if (INPUT.func_b == 5)
					{
						assert(INPUT.system=="water" or INPUT.system=="hydroxide" or INPUT.system=="hydronium");
						bool in_accept_water = false;
						for (int iwater=0; iwater<water[ia].naccept; iwater++)
						{
							int water_index = water[ia].acceptO[iwater];
							if (INPUT.ele2 == "O")
							{
								if (water_index == ia2){in_accept_water = true;} 
							}
							else if (INPUT.ele2 == "H")
							{
								for (int iH=0; iH<water[water_index].nH; iH++)
								{
									if (ia2 == water[water_index].indexH[iH])
									{
										in_accept_water = true;
									}
								}
							}
						}
						if(!in_accept_water)
						{
							continue;
						}
					}

					// mohan added if the system is hydronium or hydroxide on 2019-03-14
					if(it2==ito)
					{
						bool skip = PDF_ADDED::water_ions(cel, water, io_of_ion, ia2, should_count);
						if(skip) continue;
					}
					else if(it2==itc)
					{
						should_count = true;
						//if(ia2!=INPUT.hindex-1) continue;
//						if(ia2!=1 and ia2!=7 and ia2!=13 and ia2!=19 and ia2!=25 and ia2!=31) continue;
//						if(ia2!=2 and ia2!=8 and ia2!=14 and ia2!=20 and ia2!=26 and ia2!=32) continue;
//						if(ia2!=3 and ia2!=9 and ia2!=15 and ia2!=21 and ia2!=27 and ia2!=33) continue;
//						if(ia2!=4 and ia2!=10 and ia2!=16 and ia2!=22 and ia2!=28 and ia2!=34) continue;
//						if(ia2!=5 and ia2!=11 and ia2!=17 and ia2!=23 and ia2!=29 and ia2!=35) continue;
						if(ia2!=6 and ia2!=12 and ia2!=18 and ia2!=24 and ia2!=30 and ia2!=36) continue;
					}


					// Honeycutt-Anderson decomposion of pdf
					bool is_neighbour=true;
					if(INPUT.HA_pdf==true)
					{
						is_neighbour=Honeycutt::search_neighbours(iat, iat2, HA.nadj, 
						HA.adj_index, HA.max_num_adj, INPUT.HA_nsn, INPUT.HA_nsb);
					}
					if(!is_neighbour) continue;


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

//					ofs_running << "ia2=" << ia2 << " is being seriously considered" << endl;


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
								dis = sqrt( dx*dx + dy*dy + dz*dz );
								if( dis <= rcut )
								{
									//test
//									if(dis<2.1)
//									{
//										cout << "ppppppppppppppppppppp" << endl;
//										cout << ia+1 << " " << ia2+1 << " " << dis << " ijk " << i << " " << j << " " << k << endl;
//										cout << cel.atom[it].pos[ia].x << " " << cel.atom[it].pos[ia].y << " " << cel.atom[it].pos[ia].z << endl;
//										cout << cel.atom[it2].pos[ia2].x << " " << cel.atom[it2].pos[ia2].y << " " << cel.atom[it2].pos[ia2].z << endl;
//									} 

									if(option==1)
									{
										which = int(dis / dr);
										//cout << which; int ok; cin >> ok;
										// 2 accounts fro i-j and j-i types
										// gr[which] += 2.00;
										gr[which] += 1.00;
										if (INPUT.pdf_nstd > 0 and index_pdf >= 0)
										{
											this->multi_pdf[index_pdf][which]++;
										}
									}
									// used to calculate the SSF.
									/*
									else if(option==2)
									{
										for(int ig=0; ig<this->ng; ++ig)
										{
											gr[ig]+=cos(ig*dg*dx)*2.0;
										}
									}
									*/
								} // end dis
							} // end k
						} // end j
					} // end i
				}// ia2 
			}// it2

			
			// return the "donating" O back to the "ion" O
			if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
			{
				ia = ia_save;
				should_count = true; // renxi added 20200413
			}


		}// ia
	}// it

	
	// delete arrays, classes, clean
	if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		delete[] water;
	}
	if (INPUT.system == "water" and INPUT.nPT > 0)
	{
		if (count_num > 0) should_count = true;
		else should_count = false;
	}

	//cout << "should_count = " << should_count << endl;
	if(should_count==true)
	{		
		this->count_geometry_number++;
		cout << "count_geometry_number = " << count_geometry_number << endl;
		if (INPUT.pdf_nstd > 0)
		{
			this->npdf_count[index_pdf]++;
		}
	}

	return;
}

void PDF::sort_multi_pdf(double** multi_pdf, int* npdf_count)
{
	// Normalize the multi_pdf.
	int n1=-1;
	int n2=-1;
	if (INPUT.ntype == 1)
	{
		n1 = INPUT.natom1;
		n2 = INPUT.natom1;
	}
	else if (INPUT.ntype == 2)
	{
		if (INPUT.ele1 == INPUT.id1)
		{
			n1 = INPUT.natom1;
		}
		else if (INPUT.ele1 == INPUT.id2)
		{
			n1 = INPUT.natom2;
		}
		if (INPUT.ele2 == INPUT.id1)
		{
			n2 = INPUT.natom1;
		}
		else if (INPUT.ele2 == INPUT.id2)
		{
			n2 = INPUT.natom2;
		}
	}
	if (INPUT.system == "hydroxide" or INPUT.system == "hydronium")
	{
		n1 = 1;
		if (INPUT.ele1 == "H" and INPUT.system == "hydronium")
		{
			n1 = 3;
		}
	}
	for (int ipdf = 0; ipdf < this->npdf; ipdf++)
	{
		assert(this->npdf_count[ipdf] > 0);
		for (int ir = 0; ir<int(INPUT.rcut/INPUT.dr); ir++)
		{
			double vv = 4*(pow(ir+1, 3) - pow(ir, 3))*PI*pow(INPUT.dr, 3)/3;
			double rho = n2/INPUT.celldm1/INPUT.celldm2/INPUT.celldm3;
			this->multi_pdf[ipdf][ir] /= this->npdf_count[ipdf]*vv*n1*rho;
			//cout << ir << " factor = " << 1/(vv*n1*rho) << endl;
		}
	}
	return;
}

void PDF::out_multiple_pdf(double** multi_pdf)
{
	for (int ipdf = 0; ipdf<this->npdf; ipdf++)
	{
		ofstream ofs_multi("pdf" + to_string(ipdf) + ".txt");
		ofs_multi << "r g(r)" << endl;
		for (int ir=0; ir<int(INPUT.rcut/INPUT.dr); ir++)
		{
			ofs_multi << (ir+0.5)*INPUT.dr << " " << multi_pdf[ipdf][ir] << endl;
		}
		ofs_multi.close();
	}
	return;
}

bool PDF::correct_ion_correct_time(const Cell &cel, Water* water, int &ito, int io_of_ion)
{
	assert(INPUT.lower_time > 0);
	assert(INPUT.upper_time > 0);
	for (int is=0; is<INPUT.nPT-1; is++)
	{
		if (INPUT.non_return > 0 and this->return_jump[is])
		{
			continue;
		}
		if (cel.snapshot_time > this->snapshot_time_pt[is] - INPUT.lower_time and cel.snapshot_time < this->snapshot_time_pt[is] - INPUT.upper_time)
		{
			if (INPUT.system == "hydroxide" or INPUT.system == "hydronium")
			{
				if (io_of_ion == this->iindex_p[is])
				{
					return true;
				}
			}
			if (INPUT.system == "water")
			{
				if (io_of_ion == this->iindex[is])
				{
					return true;
				}
			}
		}
	}
	return false;
}