#include "cellFile.h"
#include "input.h"

bool CellFile::file_open=false;
ifstream CellFile::ifs_pos_kept;
ifstream CellFile::ifs_cel_kept;
ifstream CellFile::ifs_wan_kept;
ifstream CellFile::ifs_eig_kept;
ifstream CellFile::ifs_pos_ili_kept;
ifstream CellFile::ifs_vel_kept;
bool CellFile::first_read  = true;
double CellFile:: celldm = 0.0;

CellFile::CellFile()
{
}

CellFile::~CellFile()
{
}

bool CellFile::CheckGeometry( Cell &cel )
{
	TITLE("CellFile","CheckGeometry");

	// (1) check geometry file from PROFESS
	if(INPUT.geo_in_type=="PROFESS")
	{
		return CheckGeometry_PROFESS(cel);
	}
	// (2) check geometry file from VASP
	else if(INPUT.geo_in_type=="VASP")
	{
		return CheckGeometry_VASP(cel);
	}
	// (3) check geometry file from Quantum espresso.
	else if(INPUT.geo_in_type=="QE")
	{
//		CellQE::CheckGeometry(cel_in3);
	}
	// (4) check geometry file from XYZ
	else if(INPUT.geo_in_type=="XYZ")
	{
	}
	// (5) check geometry file from LAMMPS
	else if(INPUT.geo_in_type=="LAMMPS")
	{
	//	return CheckGeometry_LAMMPS(cel);
	}
	// (6) check geometry file from ABACUS
	else if(INPUT.geo_in_type=="ABACUS")
	{
		return CheckGeometry_ABACUS(cel);
	}
	// (7) check geometry file from RAW (DeePMD-kit), mohan added 2019-03-17
	else if(INPUT.geo_in_type=="RAW")
	{
		return CheckGeometry_RAW(cel); 
	}
	else
	{
		cout << " Warning! We can only accept PROFESS/VASP/QE/LAMMPS/ABACUS input geometries." << endl;
		cout << " geo_in_type = " << INPUT.geo_in_type << endl;
		exit(0);
	}
	return false;
}	




bool CellFile::ReadGeometry( Cell &cel )
{
	TITLE("CellFile","ReadGeometry");

//	cout << "ReadGeometry!" << endl;

	// (1) read in geometry file from PROFESS
	if(INPUT.geo_in_type=="PROFESS")
	{
		return ReadGeometry_PROFESS(cel);
	}
	// (2) read in geometry file from VASP
	else if(INPUT.geo_in_type=="VASP")
	{
		return ReadGeometry_VASP(cel);
	}
	else if(INPUT.geo_in_type=="VASP")
    {
        return ReadGeometry_VASP(cel);
    }
    else if(INPUT.geo_in_type=="QE2")//qianrui
    {
        if(file_open == false)
        {
            stringstream ss;
            ss << INPUT.geo_directory;
                        cout << " ReadGeometry : " << ss.str() << endl;
                        ifs_pos_kept.open(ss.str().c_str());
                        cout << " File name is " << ss.str() << endl;
                        cout << " CellFile::file_open = " << file_open << endl;
                        if(!ifs_pos_kept)
                        {
                                cout << "Could not open the file." << endl;
                                exit(0);
                        }
            file_open = true;
        }
        return ReadGeometry_QE2(cel,ifs_pos_kept);
    }
	else if(INPUT.geo_in_type=="PWmat")//qianrui
    {
        if(file_open == false)
        {
            stringstream ss;
            ss << INPUT.geo_directory;
                        cout << " ReadGeometry : " << ss.str() << endl;
                        ifs_pos_kept.open(ss.str().c_str());
                        cout << " File name is " << ss.str() << endl;
                        cout << " CellFile::file_open = " << file_open << endl;
                        if(!ifs_pos_kept)
                        {
                                cout << "Could not open the file." << endl;
                                exit(0);
                        }
            file_open = true;
        }
        return ReadGeometry_PWmat(cel,ifs_pos_kept);
    }
	// (3) read in geometry file from Quantum espresso.
	else if(INPUT.geo_in_type=="QE")
	{
		// (1) open the .pos file, .cel file, and .wan file.
		if(file_open == false)
		{
			// .pos file
			stringstream ss;
			ss << INPUT.geo_directory;
			cout << " ReadGeometry : " << ss.str() << endl;
			ifs_pos_kept.open(ss.str().c_str());
			if(!ifs_pos_kept)
			{
				cout << "could not find the .pos file: " << INPUT.geo_directory << endl;
				exit(0);
			} 
			cout << " File name is " << ss.str() << endl;
			cout << " CellFile::file_open = " << file_open << endl;


			// .cel file
			if(INPUT.cell_file!="none")
			{
				stringstream ss2;
				ss2 << INPUT.cell_file;
				cout << "Cell File: " << ss2.str() << endl;
				ifs_cel_kept.open(ss2.str().c_str());
				if(!ifs_cel_kept)
				{
					cout << "could not find the .cel file: " << INPUT.cell_file << endl;
				}
			}

			// .wan file
			if(INPUT.wannier_file!="none")
			{
				stringstream ss3;
				ss3 << INPUT.wannier_file;
				cout << " Wannier Center File: " << ss3.str() << endl;
				ifs_wan_kept.open(ss3.str().c_str());
				if(!ifs_wan_kept)
				{
					cout << "could not find the .wfc file: " << INPUT.geo_directory << endl;
					exit(0);
				}
			}

			// .eig file
			if(INPUT.eig_file!="none")
			{
				stringstream ss4;
				ss4 << INPUT.eig_file;
				cout << " Eigenvalue File: " << ss4.str() << endl;
				ifs_eig_kept.open(ss4.str().c_str());
				if(!ifs_eig_kept)
				{
					cout << "could not find the .eig file: " << INPUT.eig_file << endl;
					exit(0);
				}
			}

			// .pos_ili file
			if(INPUT.pos_ili_file!="none")
			{
				stringstream ss5;
				ss5 << INPUT.pos_ili_file;
				cout << " Pos_ILI File: " << ss5.str() << endl;
				ifs_pos_ili_kept.open(ss5.str().c_str());
				if(!ifs_pos_ili_kept)
				{
					cout << "could not find the pos_ili.dat file: " << INPUT.pos_ili_file << endl;
					exit(0);
				}	
			}
			if(INPUT.vel_file!="none")
            {
				stringstream ss6;
                ss6 << INPUT.vel_file;
                cout << " Velocity File: " << ss6.str() << endl;
                ifs_vel_kept.open(ss6.str().c_str());
                if(!ifs_vel_kept)
                {
                    cout << "could not find the vel.dat file: " << INPUT.vel_file << endl;
                    exit(0);
                }
            }

			file_open = true;
		}
			
		return ReadGeometry_QE(cel, ifs_pos_kept, ifs_cel_kept, ifs_wan_kept, ifs_eig_kept, ifs_pos_ili_kept, ifs_vel_kept);
	}
	// (4) read in geometry file from XYZ
	else if(INPUT.geo_in_type=="XYZ")
	{
		// (1) open the .pos file and .wan file.
		if(file_open == false)
		{
			stringstream ss;
			ss << INPUT.geo_directory;
			cout << " ReadGeometry : " << ss.str() << endl;
			ifs_pos_kept.open(ss.str().c_str());
			if(!ifs_pos_kept)
			{
				cout << "could not find the .pos file: " << INPUT.geo_directory << endl;
				exit(0);
			} 
			cout << " File name is " << ss.str() << endl;
			cout << " CellFile::file_open = " << file_open << endl;
			file_open = true;
		}
		return ReadGeometry_XYZ(cel, ifs_pos_kept);
	}
	// (5) read in geometry file from LAMMPS.
	else if(INPUT.geo_in_type=="LAMMPS")
	{
		if(file_open == false)
		{
			stringstream ss;
			ss << INPUT.geo_directory;
			cout << " ReadGeometry : " << ss.str() << endl;
			ifs_pos_kept.open(ss.str().c_str());
			cout << " File name is " << ss.str() << endl;
			cout << " CellFile::file_open = " << file_open << endl;

			if(!ifs_pos_kept)
			{
				cout << "Could not open the file." << endl;
				exit(0);
			} 
			file_open = true;
		}
					
		return ReadGeometry_LAMMPS(cel, ifs_pos_kept);
	}
	// (6) read in geometry file from ABACUS
	/*
	else if(INPUT.geo_in_type=="ABACUS")
	{
		return ReadGeometry_ABACUS(cel);
	}
	*/
	else if (INPUT.geo_in_type == "ABACUS")
	{
		if (file_open == false)
		{
			stringstream ss;
			ss << INPUT.geo_directory;
			cout << " ReadGeometry : " << ss.str() << endl;
			ifs_pos_kept.open(ss.str().c_str());
			if(!ifs_pos_kept)
			{
				cout << "could not find the .pos file: " << INPUT.geo_directory << endl;
				exit(0);
			} 
			cout << " File name is " << ss.str() << endl;
			cout << " CellFile::file_open = " << file_open << endl;
			file_open = true;
		}
		return ReadGeometry_ABACUS(cel, ifs_pos_kept);
	}
	// (7) read in geometry file from RAW
	else if(INPUT.geo_in_type=="RAW")
	{
        if(file_open == false)
        {
            stringstream ss;
            ss << INPUT.geo_directory;
            cout << " ReadGeometry : " << ss.str() << endl;
            ifs_pos_kept.open(ss.str().c_str());
            cout << " File name is " << ss.str() << endl;
            cout << " CellFile::file_open = " << file_open << endl;
            if(!ifs_pos_kept)
            {
                cout << "Could not open the file." << endl;
                exit(0);
            }
            file_open = true;
        }
		return ReadGeometry_RAW(cel, ifs_pos_kept);
	}
	else
	{
		cout << " Warning! We can only accept PROFESS/VASP/QE/LAMMPS/ABACUS input geometries." << endl;
		cout << " geo_in_type = " << INPUT.geo_in_type << endl;
		exit(0);
	}
	return false;
}	

void CellFile::WriteGeometry( Cell &cel, bool cartesian )
{
	ofstream ofs( INPUT.geo_out.c_str() );

	if(INPUT.geo_out_type=="PROFESS")
	{
		cout << " write geometry for PROFESS " << endl;
		WriteGeometry_PROFESS( cel, cartesian );
	}
	else if(INPUT.geo_out_type=="VASP")
	{
		cout << " write geometry for vasp " << endl;
        WriteGeometry_VASP( cel, cartesian );
		WriteGeometry_XYZ( cel );
	}
	else if(INPUT.geo_out_type=="LAMMPS")
	{
		WriteGeometry_LAMMPS( cel );
	}
	// output geometry file in the format of RAW
	else if(INPUT.geo_out_type=="RAW")
	{
		WriteGeometry_RAW( cel );
	}
	else
	{
		cout << " Warning! We can only write PROFESS/VASP/QE/LAMMPS/ABACUS input geometries." << endl;
		cout << " geo_out_type = " << INPUT.geo_out_type << endl;
		exit(0);
	}
	return;
}


bool CellFile::ReadVelocity( Cell &cel )
{
	TITLE("CellFile","ReadVelocity");

	// (1) read in velocity file from PROFESS
	if(INPUT.velcor_in_type=="PROFESS")
	{
		return ReadVelocity_PROFESS(cel);
	}
	// (2) read in velocity file from VASP
	else if(INPUT.velcor_in_type=="VASP")
	{
	}
	// (3) read in velocity file from Quantum espresso.
	else if(INPUT.velcor_in_type=="QE")
	{
	}
	// (4) read in velocity file from LAMMPS.
	//qianrui add 2020-5-10
	else if(INPUT.velcor_in_type=="LAMMPS")
	{
		//INPUT.read_velocity=true;
		INPUT.geo_in_type="LAMMPS";
		return ReadGeometry(cel);//Read velocity while reading position.
	}
	else
	{
		cout << " Error in CellFile::ReadVelocity, please indiate velcor_in_type" << endl;
		exit(0);
	}
	return false;
}	


