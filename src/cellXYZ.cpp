#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_XYZ( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_XYZ");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_XYZ( Cell &cel, ifstream &ifs)
{
	TITLE("CellFile","ReadGeometry_XYZ");
	const int ntype = INPUT.ntype;
	bool restart = true;

	assert(INPUT.celldm1>0.0);
	assert(INPUT.celldm2>0.0);
	assert(INPUT.celldm3>0.0);

	cel.a1.x = INPUT.celldm1;
	cel.a1.y = 0.0;
	cel.a1.z = 0.0;

	cel.a2.x = 0.0;
	cel.a2.y = INPUT.celldm2;
	cel.a2.z = 0.0;

	cel.a3.x = 0.0;
	cel.a3.y = 0.0;
	cel.a3.z = INPUT.celldm3;


	// (3) calculate the volume of the cell.
	cel.cal_volume();

	//cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;

	cel.atom[0].na = INPUT.natom1;
	cel.atom[0].id = INPUT.id1;
	if(INPUT.ntype>=2) {cel.atom[1].na = INPUT.natom2; cel.atom[1].id = INPUT.id2; }
	if(INPUT.ntype>=3) {cel.atom[2].na = INPUT.natom3; cel.atom[2].id = INPUT.id3; }


	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		cel.nat+=cel.atom[it].na;
		//cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}
	cout << " Total atom number is " << cel.nat << endl;

	//-----------------------
	// begin reading atoms
	//-----------------------
	assert(INPUT.natom>0);
	int tmp_natom=0;
	ifs >> tmp_natom;
	cout << tmp_natom << endl;
	assert(tmp_natom==INPUT.natom);
	string title;
	READ_VALUE(ifs, title);
	cout << title << endl;

	int *count_atom = new int[INPUT.ntype];
	for(int it=0; it<INPUT.ntype; ++it)
	{
		count_atom[it] = 0;
	}

	for(int iat=0; iat<INPUT.natom; ++iat)
	{
		string tid;
		double tx, ty, tz;
		ifs >> tid >> tx >> ty >> tz;
	//	cout << setw(5) << tid << setw(10) << cel.atom[0].id << endl;
		for(int it=0; it<ntype; ++it)
		{
			if(tid == cel.atom[it].id)
			{
				int iaa = count_atom[it];
				cel.atom[it].pos[iaa].x = tx;
				cel.atom[it].pos[iaa].y = ty;
				cel.atom[it].pos[iaa].z = tz;

				if(INPUT.length_unit=="bohr")
				{
					cel.atom[it].pos[iaa].x *= 0.529177;
					cel.atom[it].pos[iaa].y *= 0.529177;
					cel.atom[it].pos[iaa].z *= 0.529177;
				}

				while(cel.atom[it].pos[iaa].x<0) cel.atom[it].pos[iaa].x += cel.a1.norm();
				while(cel.atom[it].pos[iaa].y<0) cel.atom[it].pos[iaa].y += cel.a2.norm();
				while(cel.atom[it].pos[iaa].z<0) cel.atom[it].pos[iaa].z += cel.a3.norm();

				while(cel.atom[it].pos[iaa].x>=cel.a1.norm()) cel.atom[it].pos[iaa].x -= cel.a1.norm();
				while(cel.atom[it].pos[iaa].y>=cel.a2.norm()) cel.atom[it].pos[iaa].y -= cel.a2.norm();
				while(cel.atom[it].pos[iaa].z>=cel.a3.norm()) cel.atom[it].pos[iaa].z -= cel.a3.norm();
				++count_atom[it];
				break;
			}
		}
	}
	for(int it=0; it<ntype; ++it)
	{
    //    cout << "count_atom=" << count_atom[it] << " na=" << cel.atom[it].na << endl;
		assert(count_atom[it] == cel.atom[it].na);
	}
	delete[] count_atom;
	// ------------------------
	// finish reading atoms
	// ------------------------

	for(int it=0; it<ntype; ++it)
	{
		for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
		{
			cel.cartesian2direct(it, ia2);
//			cout << cel.atom[it].pos[ia2].x
//				<< " " << cel.atom[it].pos[ia2].y
//				<< " " << cel.atom[it].pos[ia2].z << endl;
		}
	}

	return true;
}

void CellFile::WriteGeometry_XYZ( Cell &cel )
{
    TITLE("CellFile","WriteGeometry_XYZ");
    ofstream ofs("FinalGeometry.xyz");

    ofs << cel.nat << endl;
    ofs << cel.system_name << endl;
    for(int it=0; it<INPUT.ntype; ++it)
    {
        cout << " type=" << it+1 << endl;
        cout << " natom=" << cel.atom[it].na << endl;
        for(int ia=0; ia<cel.atom[it].na; ++ia)
        {
            ofs << cel.atom[it].id
            << " " << cel.atom[it].pos[ia].x
            << " " << cel.atom[it].pos[ia].y
            << " " << cel.atom[it].pos[ia].z << endl;
        }
    }

    ofs.close();
    return;
}
