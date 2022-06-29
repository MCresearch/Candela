#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_PROFESS( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_PROFESS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << "ion." << cel.file_name << ".dat";
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_PROFESS( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_PROFESS");
	const int ntype = INPUT.ntype;

	// (1) open the file.
	stringstream ss;
	ss << INPUT.geo_directory;
	ss << "/ion." << cel.file_name << ".dat";
//	cout << " ReadGeometry : " << ss.str() << endl;
	ifstream ifs(ss.str().c_str());

	// mohan added 2019-03-12
	cel.snapshot_index = atoi(cel.file_name.c_str());
	cel.snapshot_time = cel.snapshot_index * INPUT.dt_snapshots; 


	if(!ifs) return false;
	// many files does not exist, so we don't print
	// every file's name.
	ofs_running << " File name is " << ss.str() << endl;


	cel.atom = new Atoms[ntype];

	// (2) read lattice
	bool restart = true;
	if( SCAN_BEGIN(ifs, "LATTICE_CART", !restart) )
	{
		ifs >> cel.a1.x >> cel.a1.y >> cel.a1.z;
		ifs >> cel.a2.x >> cel.a2.y >> cel.a2.z;
		ifs >> cel.a3.x >> cel.a3.y >> cel.a3.z;
//		cout << " Cell: " << endl;
//		cout << " " << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
//		cout << " " << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
//		cout << " " << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

		// mohan added 2019-03-12
		INPUT.celldm1 = sqrt(cel.a1.x*cel.a1.x + cel.a1.y*cel.a1.y + cel.a1.z*cel.a1.z);
		INPUT.celldm2 = sqrt(cel.a2.x*cel.a2.x + cel.a2.y*cel.a2.y + cel.a2.z*cel.a2.z); 
		INPUT.celldm3 = sqrt(cel.a3.x*cel.a3.x + cel.a3.y*cel.a3.y + cel.a3.z*cel.a3.z);
	}

	// (3) calculate the volume of the cell.
	cel.volume = cel.a1.x*cel.a2.y*cel.a3.z + cel.a1.y*cel.a2.z*cel.a3.x + cel.a1.z*cel.a2.x*cel.a3.y -
	  cel.a1.x*cel.a2.z*cel.a3.y - cel.a1.y*cel.a2.x*cel.a3.z - cel.a1.z*cel.a2.y*cel.a3.x;

	 //cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;


	// (4) read in atom species and the pseudopotential file.
	if( SCAN_BEGIN(ifs, "SPECIES_POT", !restart) )
	{
		for(int it=0; it<ntype; ++it)
		{
			ifs >> cel.atom[it].id >> cel.atom[it].pot_file; 
		}

		for(int it=0; it<INPUT.ntype; ++it)
		{
			if(cel.atom[it].id=="O") cel.atom[it].mass=15.9994;
			else if(cel.atom[it].id=="D") cel.atom[it].mass=2.014;
			else if(cel.atom[it].id=="H") cel.atom[it].mass=1.0079;
			else if(cel.atom[it].id=="C") cel.atom[it].mass=12.0107;
			else if(cel.atom[it].id=="Li") cel.atom[it].mass=6.941; // mohan added 2017-12-26
			else if(cel.atom[it].id=="Al") cel.atom[it].mass=26.981539; // mohan added 2019-03-12
			else
			{
				cout << "Warning! Check the atomic mass." << endl;
				exit(0);
			}
		}

	}


	// (5) read in the atom positions.
	Vector3<double> add1,add2,add3;
	if( SCAN_BEGIN(ifs, "POSITIONS_FRAC", restart) )
	{
//		string pos_type;
//		READ_VALUE(ifs, pos_type);
//		cout << " pos_type=" << pos_type << endl;
		bool frac = true;
		for(int it=0; it<ntype; ++it)
		{
			if(ntype==1) cel.atom[it].na = INPUT.natom;
			else if(ntype==2)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
			}		
			else if(ntype==3)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
				if(it==2) cel.atom[it].na = INPUT.natom3;
			}
			else if(ntype==4)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
				if(it==2) cel.atom[it].na = INPUT.natom3;
				if(it==3) cel.atom[it].na = INPUT.natom4;
			}
			else
			{
				cout << " ntype should be no more than 4." << endl;
				exit(0);
			}

			cel.atom[it].read_pos(ifs,frac);
			for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
			{
				cel.direct2cartesian(it, ia2);
/*
				cout << " direct now is " << cel.atom[it].posd[ia2].x << 
				" " << cel.atom[it].posd[ia2].y <<
				" " << cel.atom[it].posd[ia2].z << endl;
				cout << " cartesian now is " << cel.atom[it].pos[ia2].x << 
				" " << cel.atom[it].pos[ia2].y <<
				" " << cel.atom[it].pos[ia2].z << endl;
*/
			}
		
		}
	}
	else if( SCAN_BEGIN(ifs, "POSITIONS_CART", restart) )
	{
		bool frac = false;
		for(int it=0; it<ntype; ++it)
		{

			if(ntype==1) cel.atom[it].na = INPUT.natom;
			else if(ntype==2)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
				cout << " atom number is " << cel.atom[it].na << endl;
			}		
			else if(ntype==3)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
				if(it==2) cel.atom[it].na = INPUT.natom3;
			}
			else if(ntype==4)
			{
				if(it==0) cel.atom[it].na = INPUT.natom1;
				if(it==1) cel.atom[it].na = INPUT.natom2;
				if(it==2) cel.atom[it].na = INPUT.natom3;
				if(it==3) cel.atom[it].na = INPUT.natom4;
			}
			else
			{
				cout << " ntype should be no more than 4." << endl;
				exit(0);
			}


			cel.atom[it].read_pos(ifs,frac);
			for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
			{
	//			cout << " ia2=" << ia2 << endl;
				cel.cartesian2direct(it, ia2);

				// mohan added 2019-03-04
				// to make sure all atoms are within the cell
				while(cel.atom[it].posd[ia2].x<0){cel.atom[it].posd[ia2].x+=1.0;}
				while(cel.atom[it].posd[ia2].y<0){cel.atom[it].posd[ia2].y+=1.0;}
				while(cel.atom[it].posd[ia2].z<0){cel.atom[it].posd[ia2].z+=1.0;}
				while(cel.atom[it].posd[ia2].x>=1){cel.atom[it].posd[ia2].x-=1.0;}
				while(cel.atom[it].posd[ia2].y>=1){cel.atom[it].posd[ia2].y-=1.0;}
				while(cel.atom[it].posd[ia2].z>=1){cel.atom[it].posd[ia2].z-=1.0;}
				cel.direct2cartesian(it, ia2);
			}
		}
	}

	ifs.close();	

    //cout << " finish reading " << endl;

	return true;
}


void CellFile::WriteGeometry_PROFESS( Cell &cel, bool cartesian )
{
	TITLE("CellFile","WriteGeometry_PROFESS");
    // (2.1) need the number of types of elements.

	ofs_running << "geo_out:" << INPUT.geo_out << endl;

	ofstream ofs(INPUT.geo_out.c_str());

	ofs << "%BLOCK LATTICE_CART" << endl;
	ofs << setprecision(16);
	ofs << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	ofs << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	ofs << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;
	ofs << "%END BLOCK LATTICE_CART" << endl;


    if(cartesian)
	{
		ofs << "%BLOCK POSITIONS_CART" << endl;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				ofs << cel.atom[it].id 
					<< " " << cel.atom[it].pos[ia].x
					<< " " << cel.atom[it].pos[ia].y
					<< " " << cel.atom[it].pos[ia].z << endl;
			}
		}
		ofs << "%END BLOCK POSITIONS_CART" << endl;
	}
	else
	{
		ofs << "%BLOCK POSITIONS_FRAC" << endl;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << cel.atom[it].id 
					<< " " << cel.atom[it].posd[ia].x
					<< " " << cel.atom[it].posd[ia].y
					<< " " << cel.atom[it].posd[ia].z << endl;
			}
		}
		ofs << "%END BLOCK POSITIONS_FRAC" << endl;
	}


	ofs << "%BLOCK SPECIES_POT" << endl;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		ofs << cel.atom[it].id << " " << cel.atom[it].pot_file << endl;
	}
	ofs << "%END BLOCK SPECIES_POT" << endl;


	ofs << "%BLOCK ION_OPTIMIZATION" << endl;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			ofs << "1 1 1" << endl;
		}
	}
	ofs << "%EMD BLOCK ION_OPTIMIZATION" << endl;


	ofs.close();
	return;
}

bool CellFile::ReadVelocity_PROFESS( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_PROFESS");
	const int ntype = INPUT.ntype;

	// (1) open the file.
	stringstream ss;
	ss << INPUT.velcor_directory;
	ss << "/vel." << cel.file_name << ".dat";
	ifstream ifs(ss.str().c_str());
	if(!ifs)
	{
//		cout << ss.str() << " doesn't exist" << endl;
		return false; 
	}
//	cout << " File name is " << ss.str() << endl;

	cel.atom = new Atoms[ntype];

	// (2) read lattice
	bool restart = true;
	string line;
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	READ_VALUE(ifs, line);
	for(int it=0; it<INPUT.ntype; ++it)
	{
		cel.atom[it].na = INPUT.natom;
		cel.atom[it].read_vel(ifs);
	}
	return true;
}
