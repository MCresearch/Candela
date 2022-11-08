#include "cellFile.h"
#include "input.h"
#include "gfun.h"

bool CellFile::CheckGeometry_LAMMPS( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_LAMMPS");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_LAMMPS( Cell &cel, ifstream &ifs )
{
	
	TITLE("CellFile","ReadGeometry_LAMMPS");
	const int ntype = INPUT.ntype;

	//if(RANK==0) cout << "Begin reading LAMMPS input geometry" << endl;


	//----------------------------------------------------------
	// mohan added 2018-05-21
	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	// mohan updated on 2018-12-24, msd_dt is in ps!
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt; 
	// done
	//----------------------------------------------------------

	string useless;
	
	if(INPUT.geo_format==1 or INPUT.geo_format == 4)
	{
		READ_VALUE(ifs, useless); 
		READ_VALUE(ifs, useless); 
		READ_VALUE(ifs, cel.system_name);
		//cout << cel.system_name << endl;
		//cout << " Name is " << cel.system_name << endl;
		cel.nat = 0;
		READ_VALUE(ifs, cel.nat);
		assert(cel.nat == INPUT.natom);
		//int tmp_ntype = 0;
		//READ_VALUE(ifs, tmp_ntype);
		//cout << tmp_ntype << endl;
		READ_VALUE(ifs, useless);
	}
	else if(INPUT.geo_format==2)
	{
		cel.nat = 0;
		READ_VALUE(ifs, cel.nat);
		cout << "Number of atoms is " << cel.nat << endl;
		READ_VALUE(ifs, cel.nbonds);
		READ_VALUE(ifs, cel.nangles);
		cout << "Number of bonds is " << cel.nbonds << endl;
		cout << "Number of angles is " << cel.nangles << endl;
		assert(cel.nat == INPUT.natom);
		int ntype=0;
		READ_VALUE(ifs, ntype);
		READ_VALUE(ifs, cel.nbtype);
		READ_VALUE(ifs, cel.natype);
		cout << "Number of types is " << ntype << endl;
		cout << "Number of bond types is " << cel.nbtype << endl;
		cout << "Number of angle types is " << cel.natype << endl;
	}
	else if(INPUT.geo_format==3)
	{
		getline(ifs, useless);
		cout << useless << endl;
		cel.nat = 0;
		READ_VALUE(ifs, cel.nat);
		cout << "cel.nat = " << cel.nat << endl;
		int ntype = 0;
		READ_VALUE(ifs, ntype);
		cout << "ntype = " << ntype << endl;
	}

	for(int it=0; it<ntype; ++it)
	{
		//cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}
	//cout << " Total atom number is " << cel.nat << endl;



	// (2) read in cell
	double x0, x1;
	double y0, y1;
	double z0, z1;
	double xy=0.00; // mohan add 2015-06-15
	double xz=0.00;
	double yz=0.00;

	ifs >> x0;
	READ_VALUE(ifs, x1);
	ifs >> y0;
	READ_VALUE(ifs, y1);
	ifs >> z0;
	READ_VALUE(ifs, z1);

	if(INPUT.geo_format==2)
	{
		READ_VALUE(ifs, useless);
		//cout << useless << endl;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			int index; ifs >> index;
			READ_VALUE(ifs, cel.atom[it].mass);
			//cout << "mass = " << cel.atom[it].mass << endl;
		}
	}

	/*if(RANK==0)
	{
		cout << x0 << " " << x1 << endl;
		cout << y0 << " " << y1 << endl;
		cout << z0 << " " << z1 << endl;
	}*/


	// mohan add 2015-06-15
	if( INPUT.triclinic == 1 || INPUT.geo_format==3 ) 
	{
		ifs >> xy >> xz;
		READ_VALUE(ifs, yz);
	}

	cel.a1.x = x1-x0; cel.a1.y = 0.00;  cel.a1.z = 0.00;
	cel.a2.x = xy;    cel.a2.y = y1-y0; cel.a2.z = 0.00;
	cel.a3.x = xz;    cel.a3.y = yz;    cel.a3.z = z1-z0;

//	cout << " Cell: " << endl;
//	cout << " " << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
//	cout << " " << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
//	cout << " " << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;

	// mohan added 2018-05-21
	INPUT.celldm1 = cel.a1.norm();
	INPUT.celldm2 = cel.a2.norm();
	INPUT.celldm3 = cel.a3.norm();



	// (3) calculate the volume of the cell.
	cel.cal_volume();

	//cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;

	// (4) read in atom species and the pseudopotential file.
	/*
	for(int it=0; it<ntype; ++it)
	{
		stringstream ss;
		ss << "ELE" << it+1;
		cel.atom[it].id = ss.str();
	}

	string tmp_name;
	READ_VALUE(ifs, tmp_name);
	
	cout << tmp_name << endl;
	*/

	READ_VALUE(ifs, useless);
	//cout << useless << endl;

//	cout << " This is Cartesian coordinates" << endl;
	bool frac = false;
	int atom_index;
	int type_index;


	// mohan fixed 2019-04-04
	int* aifet = new int[INPUT.natom];
	for(int it=0; it<INPUT.natom; ++it)
	{
		aifet[it] = -1; // after ++ first is 0
	}
	int iat0=0;
	for(int it=0; it<ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			aifet[iat0] = ia;
			++iat0;
		}
	}		


//	int aa_index=-1;
	for(int ia=0; ia<cel.nat; ++ia)
	{	
		string useless;
		if(INPUT.geo_format==1 || INPUT.geo_format==3 || INPUT.geo_format==4)
		{
		//	ifs >> atom_index >> type_index >> useless >> useless >> useless;
			ifs >> atom_index >> type_index;
		}
		else if(INPUT.geo_format==2)
		{
			int iangle=0;
			double charge=0;
			ifs >> atom_index >> iangle >> type_index >> charge; 
			cel.atom[type_index-1].charge = charge;
			//cout << atom_index << " " << iangle << " " << type_index << " " << charge << endl;
		}
		--atom_index; // because this starts from 1 in file but should from 0 in code
		--type_index; // because this starts from 1
		//cout << type_index << endl;
		int ia0 = atom_index;
		
		
		// mohan added 2018-12-24
		if(INPUT.cartesian==true and INPUT.geo_format==1)
		{
			ifs >> cel.atom[type_index].pos[aifet[ia0]].x 
			>> cel.atom[type_index].pos[aifet[ia0]].y; 
			if(INPUT.read_velocity)//qianrui add 2020-5-20
			{
				ifs>>cel.atom[type_index].pos[aifet[ia0]].z>>cel.atom[type_index].vel[aifet[ia0]].x>>cel.atom[type_index].vel[aifet[ia0]].y;
				READ_VALUE(ifs, cel.atom[type_index].vel[aifet[ia0]].z);
			}
			else 
				READ_VALUE(ifs, cel.atom[type_index].pos[aifet[ia0]].z);
			cel.cartesian2direct(type_index, aifet[ia0]);
		}
		if (INPUT.cartesian==true and INPUT.geo_format==4)
		{
			if (type_index == 0)
			{
				ifs >> cel.atom[type_index].pos[int(ia0/3)].x 
				>> cel.atom[type_index].pos[int(ia0/3)].y;
				if(INPUT.read_velocity)//qianrui add 2020-5-20
				{
					ifs>>cel.atom[type_index].pos[int(ia0/3)].z>>cel.atom[type_index].vel[aifet[ia0]].x>>cel.atom[type_index].vel[aifet[ia0]].y;
					READ_VALUE(ifs, cel.atom[type_index].vel[aifet[ia0]].z);
				}
				else 
					READ_VALUE(ifs, cel.atom[type_index].pos[int(ia0/3)].z);
			}
			else if (type_index == 1)
			{
				int Hindex = int(ia0/3)*2;
				if (ia0%3 == 2) Hindex++;
				ifs >> cel.atom[type_index].pos[Hindex].x 
				>> cel.atom[type_index].pos[Hindex].y;
				if(INPUT.read_velocity)//qianrui add 2020-5-20
				{
					ifs>>cel.atom[type_index].pos[Hindex].z>>cel.atom[type_index].vel[aifet[ia0]].x>>cel.atom[type_index].vel[aifet[ia0]].y;
					READ_VALUE(ifs, cel.atom[type_index].vel[aifet[ia0]].z);
				}
				else 
					READ_VALUE(ifs, cel.atom[type_index].pos[Hindex].z);
			}
		}
		else if(INPUT.cartesian==false)
		{
			ifs >> cel.atom[type_index].posd[aifet[ia0]].x 
			>> cel.atom[type_index].posd[aifet[ia0]].y; 
			if(INPUT.read_velocity)
			{
				ifs>>cel.atom[type_index].posd[aifet[ia0]].z>>cel.atom[type_index].vel[aifet[ia0]].x>>cel.atom[type_index].vel[aifet[ia0]].y;
				READ_VALUE(ifs, cel.atom[type_index].vel[aifet[ia0]].z);
			}
			else 
				READ_VALUE(ifs, cel.atom[type_index].posd[aifet[ia0]].z);
			cel.direct2cartesian(type_index, aifet[ia0]);

		}

	}

	delete[] aifet;
	cel.atom_mass();

	return true;
}


void CellFile::WriteGeometry_LAMMPS( Cell &cel )
{
	TITLE("CellFile","WriteGeometry_LAMMPS");

	ofstream ofs(INPUT.geo_out.c_str());
	ofs << setprecision(16);

	ofs << "COMMENT" << endl;
	ofs << endl;

	if(INPUT.geo_format==1)
	{
		ofs << cel.nat << " atoms" << endl;
		ofs << endl;
		ofs << INPUT.ntype << " atom types" << endl;
		ofs << endl;
		ofs << "0 " << cel.a1.x << " xlo xhi" << endl;
		ofs << "0 " << cel.a2.y << " ylo yhi" << endl;
		ofs << "0 " << cel.a3.z << " zlo zhi" << endl;
		// mohan add 2015-06-15
		if( INPUT.triclinic == 1 ) 
		{
			ofs << cel.a2.x << " " << cel.a3.x << " " << cel.a3.y << " xy xz yz" << endl;
		} 
	}
	else if(INPUT.geo_format==2)
	{
		ofs << cel.nat << " atoms" << endl;
		ofs << cel.nbonds << " bonds" << endl;	
		ofs << cel.nangles << " angles" << endl;
		ofs << endl;
		ofs << INPUT.ntype << " atom types" << endl;
		ofs << "1 bond types" << endl;
		ofs << "1 angle types" << endl;
		ofs << endl;
		ofs << "0 " << cel.a1.x << " xlo xhi" << endl;
		ofs << "0 " << cel.a2.y << " ylo yhi" << endl;
		ofs << "0 " << cel.a3.z << " zlo zhi" << endl;
		ofs << endl;
		ofs << "Masses" << endl << endl;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			ofs << it+1 << " " << cel.atom[it].mass << endl;
		}
	}
	else if(INPUT.geo_format==3)
	{
		ofs << "# LAMMPS data" << endl;
		ofs << cel.nat << " atoms" << endl;
		ofs << cel.ntype << " atom types" << endl;
		ofs << "0 " << cel.a1.x << " xlo xhi" << endl;
		ofs << "0 " << cel.a2.y << " ylo yhi" << endl;
		ofs << "0 " << cel.a3.z << " zlo zhi" << endl;
		ofs << "0 0 0 xy xz yz" << endl;
		ofs << endl;
		ofs << "Atoms # metal" << endl;
		ofs << endl;
	}
	else
	{
		cout << "something wrong, chcek geo_format!" << endl;
		exit(0);
	}


	assert(cel.a1.y == 0.0);
	assert(cel.a1.z == 0.0);
	assert(cel.a2.z == 0.0);

	int ito=-1;
	int ith=-1;
	int itc=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
	}
	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
	if(INPUT.ntype==3){ assert(itc>=0); }
	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	bool frac = true;
	int iat=0;
	if(INPUT.geo_format==1)
	{
		for(int it=0; it<INPUT.ntype; ++it)
		{
			cout << " printing geometry for type " << it+1 << endl;
			cout << " there are " << cel.atom[it].na << " atoms for thie type of atom" << endl;
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				//ofs << iat+1 << " " << it+1 
				ofs << iat+1 << " " << "change_me" << it+1 // mohan update 2015-06-15 
					<< " " << cel.atom[it].pos[ia].x
					<< " " << cel.atom[it].pos[ia].y
					<< " " << cel.atom[it].pos[ia].z << endl;
				++iat;
			}
		}
	}
	else if(INPUT.geo_format==2)
	{	
		if(INPUT.system=="water")
		{
			cel.Print_water(ofs, ito, ith, norm1, norm2, norm3, INPUT.rcut_oh,1);
			cel.Print_water(ofs, ito, ith, norm1, norm2, norm3, INPUT.rcut_oh,2);
			cel.Print_water(ofs, ito, ith, norm1, norm2, norm3, INPUT.rcut_oh,3);
		}
	}
	else if(INPUT.geo_format==3)
	{
		for(int it=0; it<INPUT.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				cel.direct2cartesian(it, ia);
				ofs << iat+1 << " " << it+1 // mohan update 2015-06-15 
					<< " " << cel.atom[it].pos[ia].x
					<< " " << cel.atom[it].pos[ia].y
					<< " " << cel.atom[it].pos[ia].z << endl;
				++iat;
			}
		}
	}

	ofs.close();
	return;
}


