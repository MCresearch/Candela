#include "cellFile.h"
#include "input.h"
#include "vacuum.h"

void Vacuum::Routine()
{
	// (1) Initialization
	// cel_in : input geometry file
	// cel_out: output geometry file
	CellFile cel_in;
	CellFile cel_out;

        cel_in.file_name = INPUT.geo_in;
	
	// (2) Read in the geometry
	CellFile::ReadGeometry( cel_in );

	// (3) Extend the cell.
	Vacuum::AddVacuum( cel_in, cel_out );

	// (4) Output the geometry.
	CellFile::WriteGeometry( cel_out );

	return;
}

void Vacuum::AddVacuum( const Cell &cel_in, Cell &cel_out )
{
	// (1) initialization
	const int ntype = INPUT.ntype;
	cel_out.atom = new Atoms[ntype];

	// (2) assert we only has a1.x, a2.y and a3.z !
	assert( cel_out.a1.y == 0.0 );
	assert( cel_out.a1.z == 0.0 );
	assert( cel_out.a2.x == 0.0 );
	assert( cel_out.a2.z == 0.0 );
	assert( cel_out.a3.x == 0.0 );
	assert( cel_out.a3.y == 0.0 );

	// (3) 
	double remain_x2 = 0.0;
	double remain_y2 = 0.0;
	double remain_z2 = 0.0;
	// find the maximal value of cartesian coordinates.
	for(int it=0; it<ntype; ++it)
	{
		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{
			remain_x2 = max(remain_x2, cel_in.atom[it].pos[ia].x);
			remain_y2 = max(remain_y2, cel_in.atom[it].pos[ia].y);
			remain_z2 = max(remain_z2, cel_in.atom[it].pos[ia].z);
		}
	}
	// find the remain space in original cell.
	remain_x2 = cel_in.a1.x - remain_x2;
	remain_y2 = cel_in.a2.y - remain_y2;
	remain_z2 = cel_in.a3.z - remain_z2;

	// however, if vacuum1 and vacuum2 are both zero, means we use
	// periodic boundary condition in this direction, so, we still
	// need the original spacing between the two boundary atoms!!
	if( INPUT.vacuum_x2 == 0.0 && INPUT.vacuum_x1 == 0.0) remain_x2 = 0.0;
	if( INPUT.vacuum_y2 == 0.0 && INPUT.vacuum_y1 == 0.0) remain_y2 = 0.0;
	if( INPUT.vacuum_z2 == 0.0 && INPUT.vacuum_z1 == 0.0) remain_z2 = 0.0;
		

	// (4) Add the vacuum to the lattice vectors
	cel_out.a1.x = cel_in.a1.x + INPUT.vacuum_x1 + INPUT.vacuum_x2 - remain_x2;
	cel_out.a2.y = cel_in.a2.y + INPUT.vacuum_y1 + INPUT.vacuum_y2 - remain_y2;
	cel_out.a3.z = cel_in.a3.z + INPUT.vacuum_z1 + INPUT.vacuum_z2 - remain_z2;

	// (5) calculate the output atom positions.
	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].id = cel_in.atom[it].id;
		cel_out.atom[it].pot_file = cel_in.atom[it].pot_file;
		cel_out.atom[it].na = INPUT.natom;
		cel_out.atom[it].posd = new Vector3<double>[cel_out.atom[it].na]; 
		cel_out.atom[it].pos = new Vector3<double>[cel_out.atom[it].na]; 

		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{
			// initialize the atom positions of new cell.
			cel_out.atom[it].pos[ia].x = cel_in.atom[it].pos[ia].x + INPUT.vacuum_x1;
			cel_out.atom[it].pos[ia].y = cel_in.atom[it].pos[ia].y + INPUT.vacuum_y1; 
	          	cel_out.atom[it].pos[ia].z = cel_in.atom[it].pos[ia].z + INPUT.vacuum_z1;
			
			// change the atom positions of new cell to direct.
			cel_out.cartesian2direct(it, ia);

			//cout << cel_out.atom[it].posd[ia].x << " " << cel_out.atom[it].posd[ia].y << " " << cel_out.atom[it].posd[ia].z << endl;
		}
	}
	return;
}
