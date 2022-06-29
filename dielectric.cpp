#include "dielectric.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Dielectric::Dielectric() 
{
}

Dielectric::~Dielectric() 
{
}


void Dielectric::Routine()
{
	this->natom = INPUT.natom1; assert(natom>0);
	// compute dipole pair correlation function
	dipole_rdf();	

	return;
}

void Dielectric::dipole_rdf()
{
	TITLE("Dielectric","dipole_rdf");
	cout << "Compute dipolar pair distribution functions." << endl;

	// print out distance-dipole_angles 
    assert(INPUT.nx>0);
    assert(INPUT.ny>0);
    int nx = INPUT.nx;
    int ny = INPUT.ny;
    this->coord_xy = new double*[nx];
    for(int ix=0; ix<nx; ++ix)
    {
        this->coord_xy[ix] = new double[ny]();
    }


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
	for(int i=0; i<nmesh; ++i) gr[i] = 0.0;

	// compute Gk(R)=(1+rho*sum_0^R *cm(r)/<mu^2>)
	this->GkR = new double[nmesh];
	for(int i=0; i<nmesh; ++i) GkR[i]=0.0;

	ofs_running << " dr = " << dr << " Angstrom" << endl;
	ofs_running << " rcut = " << rcut << " Angstrom" << endl;
	ofs_running << " nmesh = " << nmesh << endl;

	// open the dipole file
	ifstream ifs(INPUT.dipole_file.c_str());
	if(!ifs)
	{
		cout << "could not find the all_dipole.dat file, quit." << endl;
		exit(0);
	}
	else
	{
		cout << "find the all_dipole.dat file: " << INPUT.dipole_file << endl;
	}
	
	// circle for all the snapshots
	assert(INPUT.geo_interval>0);
	this->count_geometry_number=0;

	cout << INPUT.geo_1 << " " << INPUT.geo_2 << endl;

	// read in the dipole moments for this snapshot.
	this->dipole_x = new double[natom]();
	this->dipole_y = new double[natom]();
	this->dipole_z = new double[natom]();
	this->dipole_tot = new double[natom]();

	// macroscopic moment
	this->M2_inter=0.0;
	this->M2_intra=0.0;
	this->M2=0.0;

	// atomic density
	this->rho_ion=0.0;
	this->volume=0.0;
	double average_dipole=0.0;

    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		CellFile cel;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) 
		{
			cout << "Error during reading geometry. Quit." << endl;
			exit(0);
		}

		int snaps=0;
		double ttt=0.0;
		double sumd=0.0;

		ifs >> snaps >> ttt;
		cout << snaps << " " << ttt << endl;
		for(int i=0; i<natom; ++i)
		{
			ifs >> dipole_x[i] >> dipole_y[i] >> dipole_z[i] >> dipole_tot[i]; 
			average_dipole += dipole_tot[i];
			sumd += dipole_tot[i];
			ofs_running << i << " " << dipole_tot[i] << endl;
		}
		ofs_running << "average dipole: " << sumd/natom << endl;
		
		this->pairs( cel, gr );

		this->rho_ion += this->natom / cel.volume;
		this->volume += cel.volume;

        cel.clean();
	}


#ifdef __MPI
	double tmp=M2_inter;
	M2_inter=0.0;
	MPI_Allreduce(&tmp, &M2_inter, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	tmp=M2_intra;
	M2_intra=0.0;
	MPI_Allreduce(&tmp, &M2_intra, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	tmp=M2;
	M2=0.0;
	MPI_Allreduce(&tmp, &M2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double *gr_local = new double[nmesh];
	for(int ir=0; ir<nmesh; ++ir)
	{
		gr_local[ir] = gr[ir];
		gr[ir] = 0.0;
	}
	MPI_Allreduce(gr_local, gr, nmesh, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	delete[] gr_local;
#endif

	// calculate the molecular number density
	//this->rho_ion = natom / (12.44*12.44*12.44);

	const double prec = 4.0/3.0*PI;
//	double unit = 0.20819434*1.602177*1.0e-19*1.0e-10; // in C*m
    if(count_geometry_number>0)
    {
		this->rho_ion /= (double)count_geometry_number;
		cout << "average ion density is " << rho_ion << " number of atoms per Angstroms^3" << endl;
		ofs_running << "average ion density is " << rho_ion << " number of atoms per Angstroms^3" << endl;
        for(int ir=0; ir<nmesh; ++ir)
        {
            gr[ir] /= (count_geometry_number*natom);
			// volume in A^3
			double vv = prec*(pow((ir+1)*dr,3)-pow(ir*dr,3));
			// rho_ion in A^3
			gr[ir] /= (rho_ion*vv);
//			gr[ir] = gr[ir]*unit*unit/EPSILON0/4/PI;
        }
		this->M2_inter /= (double)count_geometry_number;
		this->M2_intra /= (double)count_geometry_number;
		this->M2 /= (double)count_geometry_number;

		if(RANK==0)
		{
			cout << "average M2= " << M2 << " Debye^2" << endl;
			ofs_running << "average M2= " << M2 << " Debye^2" << endl;
			ofs_running << "M2_intra= " << M2_intra << endl;
			ofs_running << "M2_inter= " << M2_inter << endl;
		}

		this->volume /= (double)count_geometry_number;
		average_dipole /= ((double)count_geometry_number*natom);
    }



	// output the final pair distribution function
	if(RANK==0)
	{
		cout << "average_dipole is " << average_dipole << " Debye" << endl; 
		ofs_running << "average_dipole is " << average_dipole << " Debye" << endl; 
		cout << "the length unit of dipole_rdf.dat is Bohr" << endl;
		ofstream ofs("dipole_rdf.dat");
		assert(INPUT.mu!=0.0);
		// we output the pair correlation function and static structrue factor
		for(int i=0; i<nmesh-1; ++i)
		{
			ofs << (i+0.5)*dr/BOHR << " " << gr[i]/INPUT.mu/INPUT.mu  << endl;
		}
		ofs.close();
	}

	if(RANK==0)
	{
		this->compute_GkR(gr);
	}


	// coord_xy
	ofstream ofs2D("dis_dipole_angle_2D.dat");

    for(int ix=0; ix<nx; ++ix)
    {
		double sum=0.0;
    	for(int iy=0; iy<ny; ++iy)
        {
			sum += coord_xy[ix][iy] * INPUT.dy;
		}
	
		if(sum>0.0)
		{
			for(int iy=0; iy<ny; ++iy)
			{
				coord_xy[ix][iy]/=sum;
			}
		}
	}

	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			ofs2D << coord_xy[ix][iy] << " ";
		}
		ofs2D << endl;
	}

    for(int ix=0; ix<nx; ++ix)
    {
        delete[] coord_xy[ix];
    }
    delete[] coord_xy;
	
	ofs2D.close();


	// clean up
	ifs.close();

	delete[] gr;
	delete[] GkR;

	delete[] dipole_x;
	delete[] dipole_y;
	delete[] dipole_z;
	delete[] dipole_tot;

	return;
}

void Dielectric::pairs( const Cell &cel, double *gr)
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
	//ofs_running << "ncell_1,2,3: " << ncell_1 << " " << ncell_2 << " " << ncell_3 << endl;

	int ito=-1;
	int ith=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
	}

	bool within_distance;
	double dx0, dy0, dz0; // difference of coordinates between atom 1 and atom2
	double dx, dy, dz; // delta x,y,z between atom1, atom2.
	double x2, y2, z2;
	int which=0;
	double M2_tmp=0.0;

	// intra molecular contributions
	double M2_intra0=0.0;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		M2_intra0 += dipole_x[ia]*dipole_x[ia]+dipole_y[ia]*dipole_y[ia]+dipole_z[ia]*dipole_z[ia];
	}
	this->M2_intra += M2_intra0;

	// inter molecular contributions
	double M2_inter0=0.0;

#ifdef __MPI
	int num_atom_init = 0;
	int num_atom = cel.atom[ito].na/NPROC;
	int num_atom_extra = cel.atom[ito].na%NPROC;
	if(RANK<num_atom_extra) 
	{
		num_atom++;
		num_atom_init = RANK*num_atom;
	}
	else if(RANK>=cel.atom[ito].na)
	{
		num_atom=0;
	}
	else
	{
		num_atom_init = num_atom_extra + RANK*num_atom;
	}

	ofs_running << " num_atom_init: " << num_atom_init << endl;
	ofs_running << " num_atom: " << num_atom << endl;
	ofs_running << " num_atom_extra: " << num_atom_extra << endl;

	for(int ia=num_atom_init; ia<num_atom_init+num_atom; ++ia)
#else
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
#endif
	{
		for(int ia2=ia+1; ia2<cel.atom[ito].na; ++ia2)
		{

			// check x coordinates
			within_distance = false;
			dx0 = cel.atom[ito].pos[ia2].x - cel.atom[ito].pos[ia].x;
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
			dy0 = cel.atom[ito].pos[ia2].y - cel.atom[ito].pos[ia].y;
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
			dz0 = cel.atom[ito].pos[ia2].z - cel.atom[ito].pos[ia].z;
			for(int i=-ncell_3; i<=ncell_3; ++i)
			{
				dz = abs( dz0 + i*cel.a3.z );
				if( dz < rcut )
				{
					within_distance = true;
				}
			}
			if(within_distance == false) continue;


			for(int i=-ncell_1; i<=ncell_1; ++i)
			{
				for(int j=-ncell_2; j<=ncell_2; ++j)
				{
					for(int k=-ncell_3; k<=ncell_3; ++k)
					{
						// add cell length
						cel.add_cell_length(ito, ia2, i, j, k, x2, y2, z2);
						// calculate the distance between two atoms |r_1 - r_2|
						dx = cel.atom[ito].pos[ia].x - x2;
						dy = cel.atom[ito].pos[ia].y - y2;
						dz = cel.atom[ito].pos[ia].z - z2;
						double dis = sqrt(dx*dx+dy*dy+dz*dz);
						if(dis < rcut)
						{
							// 2 accounts fro i-j and j-i types
							which = int(dis / dr);
							M2_tmp = 2.0*(dipole_x[ia]*dipole_x[ia2]+dipole_y[ia]*dipole_y[ia2]+dipole_z[ia]*dipole_z[ia2]); 
							gr[which] += M2_tmp;	
							M2_inter0 += M2_tmp;

							// coord_xy
							int indexX = (dis-INPUT.x0)/INPUT.dx;
							//int indexY = acos(M2_tmp/2.0/dipole_tot[ia]/dipole_tot[ia2])/3.1415926*180/INPUT.dy; 
				//			cout << "cos(theta)=" << M2_tmp/2.0/dipole_tot[ia]/dipole_tot[ia2] << endl;
						
	
							int indexY = (M2_tmp/2.0/dipole_tot[ia]/dipole_tot[ia2]-INPUT.y0)/INPUT.dy; 
							// cos(theta)=d1*d2/|d1| |d2|
							if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0)
							{
								coord_xy[indexX][indexY]+=1.0;
							}
						}
					}// end k
				}// end j
			}// end i
		}// end ia2
	}// end ia

	this->count_geometry_number++;

	this->M2_inter += M2_inter0;
	
	this->M2 += M2_inter0 + M2_intra0;
	//ofs_running << "M2= " << M2 << endl;

	return;
}

void Dielectric::compute_GkR(double* gr)
{
	TITLE("Dielectric","compute_GkR");

	ofs_running << "Compute GkR based on M2_intra, rho_ion, cm(r), and dr" << endl;
	ofs_running << "M2_intra= " << M2_intra << " Dyber^2" << endl;
	ofs_running << "mu^2= " << M2_intra/this->natom << " Dyber^2" << endl;
	ofs_running << "rho_ion= " << rho_ion << " number of molecules per Angstroms^3" << endl;
	ofs_running << "dr= " << dr << " Angstrom" << endl;

	ofstream ofs("GkR.dat");

	for(int ir=0; ir<nmesh-1; ++ir)
	{
		for(int jr=0; jr<ir; ++jr)
		{
			//this->GkR[ir] += gr[jr] * this->dr * this->rho_ion;	
			this->GkR[ir] += gr[jr] * this->dr * this->rho_ion * 4 * PI; // probably wrong	
		}
//		this->GkR[ir] = this->GkR[ir]/(this->M2_intra/INPUT.natom) + 1.0;
		this->GkR[ir] = this->GkR[ir] + 1.0;

		// in unit of Bohr
		ofs << (ir+0.5)*dr/BOHR << " " << this->GkR[ir]  << endl;
	}

	ofs_running << "In order to compute the dielectric constant" << endl; 
	ofs_running << "Cell volume is " << this->volume << " Angstrom^3" << endl;

	double CM=0.20819434*1.602177*1.0e-19*1.0e-10;
	double temperature = INPUT.temperature;
	double gku=2.2;
	double M2= INPUT.mu * INPUT.mu * natom * gku * (CM * CM);
	M2/=(4*PI*EPSILON0);
	M2=1.0 + 4 * PI / 3 * M2 / (this->volume * 1.0e-30) / KB / temperature;
	ofs_running << "dielectric constant= " << M2 << endl;


	// for testing purpose only
	/*
	double sum=0.0;
	const double prec = 4.0/3.0*PI;
	for(int ir=0; ir<nmesh-1; ++ir)
	{
		double vv = prec*(pow((ir+1)*dr,3)-pow(ir*dr,3));
		sum += gr[ir]*natom*(rho_ion*vv);
	}
	ofs_running << "integral to M2_inter= " << sum << endl;

	sum=0.0;
	for(int ir=0; ir<nmesh-1; ++ir)
	{
		double rrr = ir*dr;
		sum += gr[ir]*rrr*rrr*4*PI*rho_ion;
	}
	ofs_running << "integral from rho_ion * integral_r g(r)r^2 * 4PI= " << sum << endl;
	*/

	// clean up
	ofs.close();

	return;
}
