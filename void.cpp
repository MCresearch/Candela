#include "cellFile.h"
#include "input.h"
#include "void.h"

void Void::Routine()
{
	TITLE("Void","Routine");
	// (1) Initialization
	// cel_in : input geometry file
	// cel_out: output geometry file
	CellFile cel_in;
	CellFile cel_out;

	cel_in.file_name = INPUT.geo_in; 
	
	// (2) Read in the geometry
	if( !CellFile::ReadGeometry( cel_in ) )
	{
		cout << " Can't find the geometry file : " << INPUT.geo_in << endl;
		exit(0);
	}


	// (3) Extend the cell.
	CreateVoids( cel_in, cel_out );

	bool cartesian = false;

	cartesian = INPUT.cartesian;

	// make sure this value is set.
	cel_in.ntype = INPUT.ntype;
	cel_out.ntype = cel_in.ntype;

	// (4) Output the geometry.
	cout << " should write geometry here" << endl;
	CellFile::WriteGeometry( cel_out, cartesian );

	return;
}

void Void::CreateVoids( const Cell &cel_in, Cell &cel_out )
{
	// (1) initialization
	cel_out.system_name = cel_in.system_name;
	cel_out.coordinate = cel_in.coordinate;
	cout << " cel_in.system_name=" << cel_in.system_name << endl;
	
	const int ntype = INPUT.ntype;
	cel_out.atom = new Atoms[ntype];

	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].mass = cel_in.atom[it].mass;
		cel_out.atom[it].charge = cel_in.atom[it].charge;
	}

	// (2) expand the lattice vectors
	cel_out.a1.x = cel_in.a1.x;
	cel_out.a1.y = cel_in.a1.y;
	cel_out.a1.z = cel_in.a1.z;

	cel_out.a2.x = cel_in.a2.x;
	cel_out.a2.y = cel_in.a2.y;
	cel_out.a2.z = cel_in.a2.z;

	cel_out.a3.x = cel_in.a3.x;
	cel_out.a3.y = cel_in.a3.y; 
	cel_out.a3.z = cel_in.a3.z;

	const int itref = 0;
	cout << INPUT.snatom << endl;
	assert(INPUT.snatom>0);
	double* x1 = new double[INPUT.snatom];
	double* y1 = new double[INPUT.snatom];
	double* z1 = new double[INPUT.snatom];

	int in_sphere=0;

	for(int is=0; is<INPUT.snatom; ++is)
	{
		int indexs = INPUT.satom[is];
		cout << "indexs=" << indexs << endl;
		assert(indexs>=0);
		x1[is] = cel_in.atom[itref].pos[indexs].x;
		y1[is] = cel_in.atom[itref].pos[indexs].y;
		z1[is] = cel_in.atom[itref].pos[indexs].z;
	}

	for(int it=0; it<ntype; ++it)
	{
		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{
			double x2 = cel_in.atom[it].pos[ia].x;
			double y2 = cel_in.atom[it].pos[ia].y;
			double z2 = cel_in.atom[it].pos[ia].z;
			for(int is=0; is<INPUT.snatom; ++is)
			{
				double dx=shortest(x1[is],x2,cel_in.a1.norm());
				double dy=shortest(y1[is],y2,cel_in.a2.norm());
				double dz=shortest(z1[is],z2,cel_in.a3.norm());
				double dis = sqrt(dx*dx + dy*dy + dz*dz);
				if(dis < INPUT.rcut)
				{
					++in_sphere;
					break;
				}
			}
		}
	}
	cout << "atoms within sphere: " << in_sphere << endl;

	// (3) calculate the output atom positions.
	bool frac = true;
	int ia2=0;
	cel_out.nat = 0;
	cel_out.nbonds = cel_in.nbonds; 
	cel_out.nangles = cel_in.nangles;
	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].id = cel_in.atom[it].id;
		cel_out.atom[it].pot_file = cel_in.atom[it].pot_file;

		// this sentence has problem
		cel_out.atom[it].na = cel_in.atom[it].na-in_sphere;
		cout << " number of atoms in new cell = " << cel_out.atom[it].na << endl;

		cel_out.nat += cel_out.atom[it].na;
		cel_out.atom[it].posd = new Vector3<double>[cel_out.atom[it].na]; 
		cel_out.atom[it].pos = new Vector3<double>[cel_out.atom[it].na]; 

		ia2 = 0;  // mohan fix bug 2013-06-24
		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{

			// ignore atoms in sphere
			double x2 = cel_in.atom[it].pos[ia].x;
			double y2 = cel_in.atom[it].pos[ia].y;
			double z2 = cel_in.atom[it].pos[ia].z;

			bool within_sphere=false;
			for(int is=0; is<INPUT.snatom; ++is)
			{
				double dx=shortest(x1[is],x2,cel_in.a1.norm());
				double dy=shortest(y1[is],y2,cel_in.a2.norm());
				double dz=shortest(z1[is],z2,cel_in.a3.norm());

				double dis = sqrt(dx*dx + dy*dy + dz*dz);
				if(dis < INPUT.rcut)
				{
			 		within_sphere=true;
					break;	
				}
			}
			if(within_sphere) continue;

			// calculate the new atom positions using fractional coordinates.

			cel_out.atom[it].posd[ia2].x = cel_in.atom[it].posd[ia].x;
			cel_out.atom[it].posd[ia2].y = cel_in.atom[it].posd[ia].y; 
			cel_out.atom[it].posd[ia2].z = cel_in.atom[it].posd[ia].z;

			if( cel_out.atom[it].posd[ia2].x < 0.0) cel_out.atom[it].posd[ia2].x += 1.0;
			if( cel_out.atom[it].posd[ia2].y < 0.0) cel_out.atom[it].posd[ia2].y += 1.0;
			if( cel_out.atom[it].posd[ia2].z < 0.0) cel_out.atom[it].posd[ia2].z += 1.0;

			if( cel_out.atom[it].posd[ia2].x >= 1.0) cel_out.atom[it].posd[ia2].x -= 1.0;
			if( cel_out.atom[it].posd[ia2].y >= 1.0) cel_out.atom[it].posd[ia2].y -= 1.0;
			if( cel_out.atom[it].posd[ia2].z >= 1.0) cel_out.atom[it].posd[ia2].z -= 1.0;

			cel_out.direct2cartesian(it,ia2);

			++ia2;
		}
	}
	cout << " After extension, the total atom number is " << cel_out.nat << endl;

	delete[] x1;
	delete[] y1;
	delete[] z1;
	
	return;
}
