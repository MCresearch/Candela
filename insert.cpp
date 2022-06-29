#include "cellFile.h"
#include "input.h"
#include "insert.h"
#include "random.h"

void Insert::Routine()
{
	TITLE("Insert","Routine");
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
	
	// this may not be the most suitable place
	cel_in.ntype = INPUT.ntype;

	// (3) Insert the atoms into cell. 
	Insert::InsertAtoms( cel_in, cel_out );

	bool cartesian = false;

	cartesian = INPUT.cartesian;

	// (4) Output the geometry.
	cout << " should write geometry here" << endl;
	CellFile::WriteGeometry( cel_out, cartesian );

	return;
}

void Insert::InsertAtoms( const Cell &cel_in, Cell &cel_out )
{
	TITLE("Insert","InsertAtoms");

	
	// (1) initialization
	cel_out.system_name = cel_in.system_name;
	cel_out.coordinate = cel_in.coordinate;
	cout << " cel_in.system_name=" << cel_in.system_name << endl;
	
	const int ntype = INPUT.ntype;
	cel_out.ntype = ntype+1;

	// new cell has one more element. 
	cel_out.atom = new Atoms[cel_out.ntype];

	cel_out.a1.x = cel_in.a1.x; 
	cel_out.a1.y = cel_in.a1.y;
	cel_out.a1.z = cel_in.a1.z;

	cel_out.a2.x = cel_in.a2.x;
	cel_out.a2.y = cel_in.a2.y;
	cel_out.a2.z = cel_in.a2.z;

	cel_out.a3.x = cel_in.a3.x;
	cel_out.a3.y = cel_in.a3.y;
	cel_out.a3.z = cel_in.a3.z;

	// (3) calculate the output atom positions.
	Vector3<double> add1,add2,add3;
	bool frac = true;
	int ia2=0;
	cel_out.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		cel_out.atom[it].id = cel_in.atom[it].id;
		cel_out.atom[it].pot_file = cel_in.atom[it].pot_file;

		// this sentence has problem
		cel_out.atom[it].na = cel_in.atom[it].na;

		cel_out.nat += cel_out.atom[it].na;
		cel_out.atom[it].posd = new Vector3<double>[cel_out.atom[it].na]; 
		cel_out.atom[it].pos = new Vector3<double>[cel_out.atom[it].na]; 

		ia2 = 0;  // mohan fix bug 2013-06-24
		for(int ia=0; ia<cel_in.atom[it].na; ++ia)
		{
			// calculate the new atom positions using fractional coordinates.
						
			cel_out.atom[it].posd[ia2].x = cel_in.atom[it].posd[ia].x;
			cel_out.atom[it].posd[ia2].y = cel_in.atom[it].posd[ia].y; 
			cel_out.atom[it].posd[ia2].z = cel_in.atom[it].posd[ia].z;
			cel_out.direct2cartesian(it,ia2);

			++ia2;
		}
	}

	// new inserted atoms
	cel_out.atom[ntype].id = INPUT.element_new;
	cel_out.atom[ntype].pot_file = "UNKNOWN";
	cel_out.nat += INPUT.natom_new;
	cel_out.atom[ntype].posd = new Vector3<double>[INPUT.natom_new];
	cel_out.atom[ntype].pos = new Vector3<double>[INPUT.natom_new];

	// accumulated during check the nearest distance;
	cel_out.atom[ntype].na = 0;

	// for each new atom
	for(int ia=0; ia<INPUT.natom_new; ++ia)
	{

restart:
		// direct coordinates for new atoms
		double tmpdx = Random::between0and1();
		double tmpdy = Random::between0and1();
		double tmpdz = Random::between0and1();
//		cout << " new atom " << ia+1 << " " << tmpdx << " " << tmpdy << " " << tmpdz << endl;


		// check the distances between two atoms.
		for(int it=0; it<ntype+1; ++it)
		{
			// cel_out not cel_in because we consider self distance
			for(int ia2=0; ia2<cel_out.atom[it].na; ++ia2)
			{	
				double deltax, deltay, deltaz;
				
				deltax = tmpdx - cel_out.atom[it].posd[ia2].x; 
				deltay = tmpdy - cel_out.atom[it].posd[ia2].y;
				deltaz = tmpdz - cel_out.atom[it].posd[ia2].z;

				if( abs(deltax) > abs(deltax+1.0) ) deltax = deltax+1.0;
				if( abs(deltax) > abs(deltax-1.0) ) deltax = deltax-1.0;

				if( abs(deltay) > abs(deltay+1.0) ) deltay = deltay+1.0;
				if( abs(deltay) > abs(deltay-1.0) ) deltay = deltay-1.0;

				if( abs(deltaz) > abs(deltaz+1.0) ) deltaz = deltaz+1.0;
				if( abs(deltaz) > abs(deltaz-1.0) ) deltaz = deltaz-1.0;

				//cout << " delta = " << deltax << " " << deltay << " " << deltaz << endl;
				
				// distance between two atoms
				double tmpcx = deltax * cel_in.a1.x + deltay * cel_in.a2.x + deltaz * cel_in.a3.x;
				double tmpcy = deltax * cel_in.a1.y + deltay * cel_in.a2.y + deltaz * cel_in.a3.y;
				double tmpcz = deltax * cel_in.a1.z + deltay * cel_in.a2.z + deltaz * cel_in.a3.z;

				double distance = tmpcx*tmpcx + tmpcy*tmpcy + tmpcz * tmpcz;
				distance = sqrt(distance);

				//cout << " distance = " << distance << endl;
				
				if(distance < INPUT.min_dis) 
				{
					cout << " !!!!!!!!!!!!!!!!!!!!" << endl;
					cout << " atom " << ia+1 << " in elment " << INPUT.element_new << endl;
					cout << " atom " << ia2+1 << " in species " << it+1 << endl;
					cout << " the random coordinate results distance (Angstrom) " << distance << endl;
					cout << " smaller than allowed distance " << INPUT.min_dis << endl;
					goto restart;
				}
			}	
		}
		


		cel_out.atom[ntype].na += 1;
//		cout << " already find new atom " << cel_out.atom[ntype].na << endl;

		cel_out.atom[ntype].posd[ia].x = tmpdx;
		cel_out.atom[ntype].posd[ia].y = tmpdy; 
		cel_out.atom[ntype].posd[ia].z = tmpdz; 
	}

	cout << " The new element name is " << INPUT.element_new << endl;

	cout << " After Insertion, the total atom number is " << cel_out.nat << endl;
	return;
}
