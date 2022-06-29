#include "cellFile.h"
#include "input.h"
#include "ext.h"

void Extend::Routine()
{
	TITLE("Extend","Routine");

	cout << "Extend the input cell." << endl;
	ofs_running << "Extend the input cell." << endl;
	assert(INPUT.ext_1>0);
	assert(INPUT.ext_2>0);
	assert(INPUT.ext_3>0);

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
	cout << " Read geometry type done." << endl;


	// (3) Extend the cell.
	Extend::ExtendCell( cel_in, cel_out );

	bool cartesian = false;

	cartesian = INPUT.cartesian;
	ofs_running << "cartesian: " << cartesian << endl;

	// make sure this value is set.
	cel_in.ntype = INPUT.ntype;
	cel_out.ntype = cel_in.ntype;

	// (4) Output the geometry.
	cout << " Save geometry here: " << INPUT.geo_out << endl;
	CellFile::WriteGeometry( cel_out, cartesian );
	cout << " Finished." << endl;

	return;
}

void Extend::ExtendCell( const Cell &cel_in, Cell &cel_out )
{
	TITLE("Extend","ExtendCell");

	
	// (1) initialization
	cel_out.system_name = cel_in.system_name;
	cel_out.coordinate = cel_in.coordinate;
	cout << " cel_in.system_name=" << cel_in.system_name << endl;
	cout << " cell lengths: " << cel_in.a1.x << " " << cel_in.a2.y << " " << cel_in.a3.z << endl;
	
	const int ntype = INPUT.ntype;
	cel_out.atom = new Atoms[ntype];

	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].mass = cel_in.atom[it].mass;
		cel_out.atom[it].charge = cel_in.atom[it].charge;
	}

	// (2) expand the lattice vectors
	cel_out.a1.x = cel_in.a1.x * INPUT.ext_1;
	cel_out.a1.y = cel_in.a1.y * INPUT.ext_1;
	cel_out.a1.z = cel_in.a1.z * INPUT.ext_1;

	cel_out.a2.x = cel_in.a2.x * INPUT.ext_2;
	cel_out.a2.y = cel_in.a2.y * INPUT.ext_2;
	cel_out.a2.z = cel_in.a2.z * INPUT.ext_2;

	cel_out.a3.x = cel_in.a3.x * INPUT.ext_3;
	cel_out.a3.y = cel_in.a3.y * INPUT.ext_3;
	cel_out.a3.z = cel_in.a3.z * INPUT.ext_3;

	// (3) calculate the output atom positions.
	Vector3<double> add1,add2,add3;
	bool frac = true;
	int ia2=0;
	cel_out.nat = 0;
	cel_out.nbonds = cel_in.nbonds * INPUT.ext_1 * INPUT.ext_2 * INPUT.ext_3; 
	cel_out.nangles = cel_in.nangles * INPUT.ext_1 * INPUT.ext_2 * INPUT.ext_3; 
	cout << " Extension = " << INPUT.ext_1 << " * " << INPUT.ext_2 << " * " << INPUT.ext_3 << endl; 
	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].id = cel_in.atom[it].id;
		cel_out.atom[it].pot_file = cel_in.atom[it].pot_file;

		// this sentence has problem
		cel_out.atom[it].na = cel_in.atom[it].na * INPUT.ext_1 * INPUT.ext_2 * INPUT.ext_3;
		cout << " number of atoms in new cell = " << cel_out.atom[it].na << endl;

		cel_out.nat += cel_out.atom[it].na;
		cel_out.atom[it].posd = new Vector3<double>[cel_out.atom[it].na]; 
		cel_out.atom[it].pos = new Vector3<double>[cel_out.atom[it].na]; 

		ia2 = 0;  // mohan fix bug 2013-06-24
		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{
		//	cout << cel_out.atom[it].posd[ia].x << " " << cel_out.atom[it].posd[ia].y << " " << cel_out.atom[it].posd[ia].z << endl;
			for(int i=0; i<INPUT.ext_1; i++)
			{
				for(int j=0; j<INPUT.ext_2; j++)
				{
					for(int k=0; k<INPUT.ext_3; k++)
					{
						// calculate the new atom positions using fractional coordinates.
						
						cel_out.atom[it].posd[ia2].x = cel_in.atom[it].posd[ia].x + (double)i;
						cel_out.atom[it].posd[ia2].y = cel_in.atom[it].posd[ia].y + (double)j; 
						cel_out.atom[it].posd[ia2].z = cel_in.atom[it].posd[ia].z + (double)k;


						cel_out.atom[it].posd[ia2].x /= (double)INPUT.ext_1;
						cel_out.atom[it].posd[ia2].y /= (double)INPUT.ext_2;
						cel_out.atom[it].posd[ia2].z /= (double)INPUT.ext_3;

                        // mohan added 2016-05-03
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
			}
		}
	}
	cout << " After extension, the total atom number is " << cel_out.nat << endl;
	return;
}
