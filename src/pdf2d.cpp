#include "cellFile.h"
#include "input.h"
#include "pdf2d.h"
#include "math.h"

void PDF2d::Routine()
{
	TITLE("PDF2d","Routine");

	cal();

	return;
}


void PDF2d::cal()
{
	TITLE("PDF2d","cal");

	// (1) set the basic parameters.
	// delta r in real space.
	this->dr = INPUT.dr;

	// radius cutoff in real space, usually choose a/2,
	// where a is the lattice constant.
	this->rcut = INPUT.rcut;

	// number of radial mesh.
	this->nmesh = int(rcut / dr) +  1;
	
	// pair distribution function.
	double* gr = new double[nmesh];
	double* gr_tmp = new double[nmesh];
	for(int i=0; i<nmesh; ++i) gr[i] = 0;
	for(int i=0; i<nmesh; ++i) gr_tmp[i] = 0;

	cout << " dr = " << dr << " Angstrom" << endl;
	cout << " rcut = " << rcut << " Angstrom" << endl;
	cout << " nmesh = " << nmesh << endl;


	// ionic density 
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;


	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
	{
		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) continue;

		this->natom_per_layer = 0;

		CellFile cel;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
		++count_geometry_number;

		// input the reference density
		assert(INPUT.rho_ref>0);
		const double rho_ion = INPUT.rho_ref;

		if(count_geometry_number==1)
		{
			cout << " Average ion density = " << rho_ion << endl;
			// Fermi vector
			double kf = pow(3*PI*PI*rho_ion,1.0/3.0);
	        cout << " Fermi vector = " << kf << endl; 
			cout << " pi=" << PI << endl;

	        assert(kf > 0.0);
		}

		// add this 2014-05-30
		for(int i=0; i<nmesh; ++i) gr_tmp[i] = 0;
		
		// (2)
		int option = 1;
		this->periodic_pairs( cel, gr_tmp, option );


		// (3) calculate the pair distribution function
		//
		// dV, equals to 4 * pi * r^2 * dr
		assert( natom_per_layer > 0);
		const double prec = 4.0/3.0*PI;
		for(int i=0; i<nmesh; ++i)
		{
			// volume
			double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));

			if(INPUT.ntype==1)
			{
				gr_tmp[i] = gr_tmp[i] / this->natom_per_layer / rho_ion / vv;
			}
			else if(INPUT.ntype==2)
			{
				int n1=0;
				int n2=0;
				for(int it=0; it<INPUT.ntype; ++it)
				{
					if(cel.atom[it].id == INPUT.ele1) n1=cel.atom[it].na;
					if(cel.atom[it].id == INPUT.ele2) n2=cel.atom[it].na;
				}
				assert(n1>0);
				assert(n2>0);
				gr_tmp[i] = gr_tmp[i] * this->natom_per_layer / n1 /n2 / rho_ion / vv;
			}
		}

		for(int i=0; i<nmesh; ++i)
		{
			gr[i]+=gr_tmp[i];
		}
	}

    // do average of 3D ssf.
    assert(count_geometry_number>0);
    cout << " count_geometry_number = " << count_geometry_number << endl;
    if(count_geometry_number>0)
    {
        for(int ir=0; ir<nmesh; ++ir)
        {
            gr[ir] /= count_geometry_number;
        }
    }



	// output the final pair distribution function
	ofstream ofs(INPUT.geo_out.c_str());
	// we output the pair correlation function and static structrue factor
	for(int i=0; i<nmesh-1; ++i)
	{
		ofs << (i+0.5)*dr << " " << gr[i]  << endl;
	}
	ofs.close();




// calculate the SSF via FFT,
/*
	if(ng>0 && dg>0)
	{
		// calculate the static structure factor 
		this->dg = INPUT.struf_dg;
		this->ng = INPUT.struf_ng;

		double *sf = new double[ng];
		for(int ig=0; ig<ng; ++ig) sf[ig] = 0.0;

		this->static_structure_factor( cel, rho_ion, gr, sf );

		// output the static structure factor.
		ofstream ofss(INPUT.ssf_out.c_str());
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
*/
	
	// --> clean <--

	delete[] gr;
	return;
}


void PDF2d::periodic_pairs( const Cell &cel, double *gr, const int &option)
{	
	// (1) calculate the norm of each lattice vectors.
	double a1 = cel.a1.norm();
	double a2 = cel.a2.norm();
	double a3 = cel.a3.norm();
//	cout << " norms of lattice vectors = " << a1 << " " << a2 << " " << a3 << endl;

	// (2) calculate how many more cells need to be convered. 
	assert(rcut>0);
	const int ncell_1 = int(rcut/a1)+1;
	const int ncell_2 = int(rcut/a2)+1;
	// for 3D:
	//const int ncell_3 = int(rcut/a3)+1;
	// for 2D:
	const int ncell_3 = 0;
//	cout << " ncell is " << ncell_1 << " " << ncell_2 << " " << ncell_3 << endl;

	// (3) calculate the distance between atoms.
	double x2, y2, z2; // atom positions for atom 2.
	double dx, dy, dz; // delta x,y,z between atom1, atom2.
	double dis;
	int which;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(cel.atom[it].id != INPUT.ele1) continue;

		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			if( cel.atom[it].pos[ia].z < INPUT.pdf_z0 ) continue;
			if( cel.atom[it].pos[ia].z > INPUT.pdf_z1 ) continue;

			// mohan add 2015-05-29
			this->natom_per_layer = natom_per_layer + 1;
//			cout << cel.atom[it].pos[ia].x << " " 
//			<< cel.atom[it].pos[ia].y << " "
//			<< cel.atom[it].pos[ia].z << endl;

			for(int it2=0; it2<INPUT.ntype; ++it2)
			{
				if(cel.atom[it2].id != INPUT.ele2) continue;

				for(int ia2=0; ia2<cel.atom[it2].na; ++ia2)
				{
					if(it==it2 && ia==ia2) continue;
			
					if( cel.atom[it2].pos[ia2].z < INPUT.pdf_z0 ) continue;
					if( cel.atom[it2].pos[ia2].z > INPUT.pdf_z1 ) continue;

					// search for 27 cells.
					bool within_distance=false;
					for(int i=-ncell_1; i<=ncell_1; ++i)
					{
						for(int j=-ncell_2; j<=ncell_2; ++j)
						{
							for(int k=-ncell_3; k<=ncell_3; ++k)
							{
								if(within_distance) continue;

								// add cell length
								cel.add_cell_length(it2, ia2, i, j, k, x2, y2, z2);

								// calculate the distance between two atoms |r_1 - r_2|
								dx = cel.atom[it].pos[ia].x - x2;
								dy = cel.atom[it].pos[ia].y - y2;
								dz = cel.atom[it].pos[ia].z - z2;
								dis = sqrt( dx*dx + dy*dy + dz*dz );
								if( dis <= rcut )
								{
									if(option==1)
									{
										which = int(dis / dr);
										//cout << which; int ok; cin >> ok;
										// 2 accounts fro i-j and j-i types
										// gr[which] += 2.00;
										gr[which] += 1.00;
									}
									within_distance=true;
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}
