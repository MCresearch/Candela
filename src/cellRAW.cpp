#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_RAW( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_RAW");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_RAW( Cell &cel, ifstream &ifs)
{
	TITLE("CellFile","ReadGeometry_RAW");
	const int ntype = INPUT.ntype;

	//----------------------------------------------------------
    // mohan added 2020-03-24
    static int count_geometry=0;
    cel.snapshot_index = count_geometry;
    ++count_geometry;
    // mohan updated on 2020-03-24, msd_dt is in ps!
    cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt;
    // done
    //----------------------------------------------------------

	assert(INPUT.celldm1>0.0);
	assert(INPUT.celldm2>0.0);
	assert(INPUT.celldm3>0.0);

	cel.a1.x = INPUT.celldm1; cel.a1.y = 0.0;  cel.a1.z = 0.0;
	cel.a2.x = 0.0; cel.a2.y = INPUT.celldm2; cel.a2.z = 0.0;
	cel.a3.x = 0.0; cel.a3.y = 0.0; cel.a3.z = INPUT.celldm3;

	// (3) calculate the volume of the cell.
	cel.cal_volume();

	//cout << " volume of the cell is " << cel.volume << " (Angstrom^3)" << endl;

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
	int *count_atom = new int[INPUT.ntype];
	for(int it=0; it<INPUT.ntype; ++it)
	{
		count_atom[it] = 0;
	}

	double tx, ty, tz;
	for(int it=0; it<ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			ifs >> tx >> ty >> tz;
			cel.atom[it].pos[ia].x = tx;
			cel.atom[it].pos[ia].y = ty;
			cel.atom[it].pos[ia].z = tz;

			while(cel.atom[it].pos[ia].x<0) cel.atom[it].pos[ia].x += cel.a1.norm();
			while(cel.atom[it].pos[ia].y<0) cel.atom[it].pos[ia].y += cel.a2.norm();
			while(cel.atom[it].pos[ia].z<0) cel.atom[it].pos[ia].z += cel.a3.norm();

			while(cel.atom[it].pos[ia].x>=cel.a1.norm()) cel.atom[it].pos[ia].x -= cel.a1.norm();
			while(cel.atom[it].pos[ia].y>=cel.a2.norm()) cel.atom[it].pos[ia].y -= cel.a2.norm();
			while(cel.atom[it].pos[ia].z>=cel.a3.norm()) cel.atom[it].pos[ia].z -= cel.a3.norm();
		}
	}

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


void CellFile::WriteGeometry_RAW( Cell &cel )
{
    TITLE("CellFile","WriteGeometry_RAW");
    ofstream ofs("coord_d310.raw");

    for(int it=0; it<INPUT.ntype; ++it)
    {
        cout << " type=" << it+1 << endl;
        cout << " natom=" << cel.atom[it].na << endl;
        for(int ia=0; ia<cel.atom[it].na; ++ia)
        {
            ofs << " " << cel.atom[it].pos[ia].x
            << " " << cel.atom[it].pos[ia].y
            << " " << cel.atom[it].pos[ia].z;
        }
    }
	ofs << endl;

    ofs.close();
    return;
}

