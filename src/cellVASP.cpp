#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_VASP( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_VASP");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_VASP( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_VASP");
	const int ntype = INPUT.ntype;

	// (1) open the file.
	stringstream ss;
	ss << INPUT.geo_directory << "/";
	ss << cel.file_name;
	cout << " ReadGeometry : " << ss.str() << endl;


	ifstream ifs(ss.str().c_str());

	if(!ifs) return false;
	// many files does not exist, so we don't print
	// every file's name.
	cout << " File name is " << ss.str() << endl;


	cel.atom = new Atoms[ntype];
	getline(ifs, cel.system_name);
	cout << " Name is " << cel.system_name << endl;

	// (2) read lattice
	bool restart = true;
	double scaling = 1.0;
	ifs >> scaling;
	cout << " scaling=" << scaling << endl;

    ifs >> cel.a1.x >> cel.a1.y >> cel.a1.z;
	ifs >> cel.a2.x >> cel.a2.y >> cel.a2.z;
	ifs >> cel.a3.x >> cel.a3.y >> cel.a3.z;

	cel.a1.x *= scaling;
	cel.a1.y *= scaling;
	cel.a1.z *= scaling;
	cel.a2.x *= scaling;
	cel.a2.y *= scaling;
	cel.a2.z *= scaling;
	cel.a3.x *= scaling;
	cel.a3.y *= scaling;
	cel.a3.z *= scaling;

	cout << " Cell: " << endl;
	cout << " " << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	cout << " " << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	cout << " " << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

	// (3) calculate the volume of the cell.
	cel.volume = cel.a1.x*cel.a2.y*cel.a3.z + cel.a1.y*cel.a2.z*cel.a3.x + cel.a1.z*cel.a2.x*cel.a3.y -
	  cel.a1.x*cel.a2.z*cel.a3.y - cel.a1.y*cel.a2.x*cel.a3.z - cel.a1.z*cel.a2.y*cel.a3.x;

	cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;


	// (4) read in atom species and the pseudopotential file.
	for(int it=0; it<ntype; ++it)
	{
		ifs >> cel.atom[it].id; 
	}

	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		ifs >> cel.atom[it].na; 
		cel.nat+=cel.atom[it].na;
		cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}
	cout << " Total atom number is " << cel.nat << endl;

    string tmp;
	string selected_dynamics;
	getline(ifs, tmp);
	getline(ifs, selected_dynamics);
	cout << " The special line : " << selected_dynamics << endl;

	cel.coordinate = "Direct";

	// (5) read in the atom positions.
	Vector3<double> add1,add2,add3;
	if( SCAN_BEGIN(ifs, "Direct", restart) )
	{
		cout << " This is Direct coordinates" << endl;
//		string pos_type;
//		READ_VALUE(ifs, pos_type);
//		cout << " pos_type=" << pos_type << endl;
		bool frac = true;
		for(int it=0; it<ntype; ++it)
		{
			cel.atom[it].read_pos_2(ifs,frac);
			for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
			{
				cel.direct2cartesian(it, ia2);
			}
		}
	}
	
	if( SCAN_BEGIN(ifs, "Cartesian", restart) )
	{
		cout << " This is Cartesian coordinates" << endl;
		bool frac = false;
		for(int it=0; it<ntype; ++it)
		{
			cel.atom[it].read_pos_2(ifs,frac);
			for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
			{
				cel.cartesian2direct(it, ia2);
			}
		}
	}
	ifs.close();	
	return true;
}


void CellFile::WriteGeometry_VASP( Cell &cel, bool cartesian )
{
	TITLE("CellFile","WriteGeometry_VASP");
	cout << " wirte out geometry for vasp." << endl;
	if( cel.ntype != INPUT.ntype )
	{
		cout << " Be careful the printted numbe of species is not equal to the input 'ntype' " << endl;
		cout << " now number of species : " << cel.ntype << endl;
		cout << " ntype from input file : " << INPUT.ntype << endl;
	}
		

    // (2.1) need the number of types of elements.
	
	ofstream ofs(INPUT.geo_out.c_str());

	ofs << "VASP_POSCAR_FORMAT" << endl;
	ofs << "1" << endl; // scaling factor
	ofs << setprecision(16);
	ofs << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
	ofs << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
	ofs << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

	// name of the element
	for(int it=0; it<cel.ntype; ++it)
	{
		if( cel.atom[it].id == "" ) 
		{
			cel.atom[it].id = "XX";
		}
		ofs << cel.atom[it].id << " ";
	}
	ofs << endl;

	// atom number
	for(int it=0; it<cel.ntype; ++it)
	{
		ofs << cel.atom[it].na << " ";
	}
	ofs << endl;

	ofs << "Selected Dynamics" << endl;

	if(!cartesian)
	{
		ofs << "Direct" << endl;
		// atom position
		for(int it=0; it<cel.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << cel.atom[it].posd[ia].x
					<< " " << cel.atom[it].posd[ia].y
					<< " " << cel.atom[it].posd[ia].z 
					<< " T T T" 
					<< endl;
			}
		}
	}
	else
	{
		ofs << "Cartesian" << endl;
		// atom position
		for(int it=0; it<cel.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << cel.atom[it].pos[ia].x
					<< " " << cel.atom[it].pos[ia].y
					<< " " << cel.atom[it].pos[ia].z 
					<< " T T T" 
					<< endl;
			}
		}
	}

	ofs.close();
	return;
}

bool CellFile::ReadVelocity_VASP( Cell &cel )
{
	TITLE("CellFile","ReadGeometry_VASP");
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
	//cout << " File name is " << ss.str() << endl;

	cel.atom = new Atoms[ntype];

	// (2) read lattice
	bool restart = true;
	string line1,line2,line3,line4;
	READ_VALUE(ifs, line1);
	READ_VALUE(ifs, line2);
	READ_VALUE(ifs, line3);
	READ_VALUE(ifs, line4);
	for(int it=0; it<INPUT.ntype; ++it)
	{
		cel.atom[it].na = INPUT.natom;
		cel.atom[it].read_vel(ifs);
	}
	return true;
}
