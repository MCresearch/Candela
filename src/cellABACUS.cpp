#include "cellFile.h"
#include "input.h"
#include "gfun.h"
#include "string.h"


bool CellFile::CheckGeometry_ABACUS( Cell &cel )
{
	/*
    TITLE("CellFile","CheckGeometry_ABACUS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
	*/
}

bool CellFile::ReadGeometry_ABACUS( Cell &cel, ifstream &ifs )
{
	TITLE("CellFile","ReadGeometry_ABACUS");
	const int ntype = INPUT.ntype;
	bool restart = true;

	string useless;
	double lat_const;
	ifs >> useless >> useless;
	ifs >> useless >> lat_const >> useless;
	assert(useless == "Angstrom");
	ifs >> useless;
	ifs >> cel.a1.x >> cel.a1.y >> cel.a1.z;
	ifs >> cel.a2.x >> cel.a2.y >> cel.a2.z;
	ifs >> cel.a3.x >> cel.a3.y >> cel.a3.z;

	cel.a1 *= lat_const;
	cel.a2 *= lat_const;
	cel.a3 *= lat_const;

	INPUT.celldm1 = cel.a1.norm();
	INPUT.celldm2 = cel.a2.norm();
	INPUT.celldm3 = cel.a3.norm();
	cel.cal_volume();

	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt; 

	char keyword[20];
	ifs >> keyword;
	if (strcmp(keyword, "VIRIAL") == 0) 
	{
		for (int i=0; i<5; i++) getline(ifs,useless);
	}
	else if (strcmp(keyword, "INDEX") == 0)
	{
		getline(ifs,useless);
	}

	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		cel.nat+=cel.atom[it].na;
		//cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}

	for(int it=0; it<ntype; ++it)
	{
		Atoms *atomit = &cel.atom[it];
		for(int i=0; i<atomit->na; ++i)
		{
			ifs >> useless >> useless;
			ifs >> atomit->pos[i].x >> atomit->pos[i].y;
			READ_VALUE(ifs, atomit->pos[i].z);
		}

		for(int ia2=0; ia2<atomit->na; ++ia2)
		{
			//atomit->pos[ia2] *= lat_const;
			// cout << atomit->pos[ia2].x << " " << atomit->pos[ia2].y << " " << atomit->pos[ia2].z << endl;
			cel.cartesian2direct(it, ia2);
		//	cout << atomit->id << " " << atomit->pos[ia2].x << " " << atomit->pos[ia2].y << " " << atomit->pos[ia2].z << endl;
//			cout << atomit->pos[ia2].x
//				<< " " << atomit->pos[ia2].y
//				<< " " << atomit->pos[ia2].z << endl;
		}
	}
	cel.atom_mass();
	return true;
}


void CellFile::WriteGeometry_ABACUS(Cell &cel)
{
	/*
	TITLE("CellFile","WriteGeometry_ABACUS");
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
	*/
}

