#include "cellFile.h"
#include "input.h"
#include "gfun.h"


bool CellFile::CheckGeometry_QE( Cell &cel )
{
    TITLE("CellFile","CheckGeometry_QE");
    const int ntype = INPUT.ntype;

    // (1) open the file.
    stringstream ss;
    ss << INPUT.geo_directory;
    ss << cel.file_name;
    ifstream ifs(ss.str().c_str());
    if(!ifs) return false;
	else return true;
}


bool CellFile::ReadGeometry_QE( Cell &cel, ifstream &ifs, ifstream &ifs_cel, ifstream &ifs_wan, ifstream &ifs_eig, ifstream &ifs_pos_ili, ifstream &ifs_vel)
{
	TITLE("CellFile","ReadGeometry_QE");
	const int ntype = INPUT.ntype;
	bool restart = true;

	delete[] cel.atom;
	cel.atom = new Atoms[ntype];

	if(INPUT.eig_file != "none")
	{
		assert(INPUT.nbands>0);
		cel.nbands = INPUT.nbands;

		delete[] cel.eig;
		cel.eig = new double[INPUT.nbands];
		cel.read_eig(ifs_eig, INPUT.nbands);
	}
	string useless; // renxi added 20191226
	// read this information first in order to do comparison
	if (INPUT.dt_snapshots <= 0)
	{
		ifs >> cel.snapshot_index >> cel.snapshot_time; //renxi commented this out; data contained in H3O-SCAN-QE has *** rather than numbers. renxi 20191226
	}
//	ifs >> cel.snapshot_index >> cel.snapshot_time;// renxi added 20191226
	if (INPUT.dt_snapshots > 0)
	{
		ifs >> cel.snapshot_index >> useless;// renxi added 20191226
		cel.snapshot_time=cel.snapshot_index*INPUT.dt_snapshots; // renxi added 20191226 fixed on 20200730
	}
//	cout << "snapshot: " << cel.snapshot_index << " " << cel.snapshot_time << endl;
	cout << "snapshot_index " << cel.snapshot_index << endl;
        cout << "snapshot_time " << cel.snapshot_time << endl;
	
	if(INPUT.wannier_file != "none")
	{
		assert(INPUT.nbands>0);
		cel.nbands = INPUT.nbands;

		delete[] cel.wan_centers;
		cel.wan_centers = new Vector3<double>[INPUT.nbands];
		cel.read_wannier_centers(ifs_wan, INPUT.nbands);
	}
	if(INPUT.cell_file=="none")
	{
		// mohan update this on Feb.23 2018 (add e12 e13 e21 e23 e31 e32)
		if(INPUT.celldm1==0.0 and INPUT.celldm2==0.0 and INPUT.celldm3==0.0)
		{
			cel.a1.x = INPUT.e11;
			cel.a1.y = INPUT.e12;
			cel.a1.z = INPUT.e13;

			cel.a2.x = INPUT.e21;
			cel.a2.y = INPUT.e22; 
			cel.a2.z = INPUT.e23;

			cel.a3.x = INPUT.e31;
			cel.a3.y = INPUT.e32;
			cel.a3.z = INPUT.e33; 

			INPUT.celldm1 = sqrt(cel.a1.x*cel.a1.x + cel.a1.y*cel.a1.y + cel.a1.z*cel.a1.z);
			INPUT.celldm2 = sqrt(cel.a2.x*cel.a2.x + cel.a2.y*cel.a2.y + cel.a2.z*cel.a2.z); 
			INPUT.celldm3 = sqrt(cel.a3.x*cel.a3.x + cel.a3.y*cel.a3.y + cel.a3.z*cel.a3.z);
		}
		else // a, b, c are orthogonal
		{
			cel.a1.x = INPUT.celldm1;
			cel.a1.y = 0.0;
			cel.a1.z = 0.0;

			cel.a2.x = 0.0;
			cel.a2.y = INPUT.celldm2;
			cel.a2.z = 0.0;

			cel.a3.x = 0.0;
			cel.a3.y = 0.0;
			cel.a3.z = INPUT.celldm3;
		}

	}
	else
	{
		int snapshot_cel=0;
		double time_cel=0.0;
		ifs_cel >> snapshot_cel >> time_cel;
		//cout << "snapshot_cel:" << snapshot_cel << endl;
		//cout << "snapshot_index: " << cel.snapshot_index << endl;
	//	cout << "time cell: " << time_cel << endl;
		assert(snapshot_cel == cel.snapshot_index);
		assert(time_cel == cel.snapshot_time);

		// mohan update according to the inputs description on QE website, 2018-05-16
		ifs_cel >> cel.a1.x >> cel.a2.x >> cel.a3.x;
		ifs_cel >> cel.a1.y >> cel.a2.y >> cel.a3.y;
		ifs_cel >> cel.a1.z >> cel.a2.z >> cel.a3.z;
		cel.a1 *= BOHR; // units: Angstrom
		cel.a2 *= BOHR; // units: Angstrom
		cel.a3 *= BOHR; // units: Angstrom
		//cout << cel.a1.x << " " << cel.a1.y << " " << cel.a1.z << endl;
		//cout << cel.a2.x << " " << cel.a2.y << " " << cel.a2.z << endl;
		//cout << cel.a3.x << " " << cel.a3.y << " " << cel.a3.z << endl;
		INPUT.celldm1 = sqrt(cel.a1.x*cel.a1.x + cel.a1.y*cel.a1.y + cel.a1.z*cel.a1.z);
		INPUT.celldm2 = sqrt(cel.a2.x*cel.a2.x + cel.a2.y*cel.a2.y + cel.a2.z*cel.a2.z); 
		INPUT.celldm3 = sqrt(cel.a3.x*cel.a3.x + cel.a3.y*cel.a3.y + cel.a3.z*cel.a3.z);
	}


	// (3) calculate the volume of the cell.
	cel.volume = cel.a1.x*cel.a2.y*cel.a3.z + cel.a1.y*cel.a2.z*cel.a3.x + cel.a1.z*cel.a2.x*cel.a3.y -
	  cel.a1.x*cel.a2.z*cel.a3.y - cel.a1.y*cel.a2.x*cel.a3.z - cel.a1.z*cel.a2.y*cel.a3.x;

	ofs_running << cel.snapshot_index << " " << cel.snapshot_time << " volume " << cel.volume << " (Angstrom^3)"
	<< " rho(64H2O)= " << 64*18*1.6605/cel.volume << endl; // temporary code

	cel.atom[0].na = INPUT.natom1;
	cel.atom[0].id = INPUT.id1;
	if(INPUT.ntype>=2) {cel.atom[1].na = INPUT.natom2; cel.atom[1].id = INPUT.id2; }
	if(INPUT.ntype>=3) {cel.atom[2].na = INPUT.natom3; cel.atom[2].id = INPUT.id3; }
//	if(INPUT.ntype>=4) {cel.atom[3].na = INPUT.natom4; cel.atom[3].id = INPUT.id4; }

	// mohan aded 2017-09-11
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

	cel.nat = 0;
	for(int it=0; it<ntype; ++it)
	{
		cel.nat+=cel.atom[it].na;
		//cout << " Element : " << cel.atom[it].id << " Atom Number : " << cel.atom[it].na << endl;
	}
	//cout << " Total atom number is " << cel.nat << endl;


	for(int it=0; it<ntype; ++it)
	{
		cel.atom[it].read_pos_4(ifs, cel.a1, cel.a2, cel.a3);
		for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
		{
			cel.cartesian2direct(it, ia2);
		//	cout << cel.atom[it].id << " " << cel.atom[it].pos[ia2].x << " " << cel.atom[it].pos[ia2].y << " " << cel.atom[it].pos[ia2].z << endl;
//			cout << cel.atom[it].pos[ia2].x
//				<< " " << cel.atom[it].pos[ia2].y
//				<< " " << cel.atom[it].pos[ia2].z << endl;
		}
	}


    // mohan added 2017-04-08
    if(INPUT.pos_ili_file != "none")
    {
        assert(INPUT.ntype>0);
        int index_ili;
        double time_ili;
        ifs_pos_ili >> index_ili >> time_ili;
    //    cout << index_ili << " " << time_ili << endl;
        for(int it=0; it<INPUT.ntype; ++it)
        {
    //        cout << "allocate_pos_ili=" << cel.atom[it].allocate_pos_ili << endl;
            if(cel.atom[it].allocate_pos_ili==false)
            {
                cel.atom[it].pos_ili = new double[cel.atom[it].na];
                cel.atom[it].allocate_pos_ili=true;
            }
            cel.read_pos_ili(ifs_pos_ili, it);
        }
    }

	if(INPUT.vel_file != "none")
	{
		int index_vel = -1;
		double time_vel = -1;

		ifs_vel >> index_vel >> time_vel;
		assert(index_vel == cel.snapshot_index);
		for(int it=0; it<INPUT.ntype; ++it)
        {
			//delete[] cel.atom[it].vel;
			cel.atom[it].read_vel(ifs_vel);
		}
	}



//	cout << "Finish reading cells." << endl;

	return true;
}


