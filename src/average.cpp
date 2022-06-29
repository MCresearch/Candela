#include "cellFile.h"
#include "input.h"
#include "average.h"
#include "math.h"

void Average::Routine()
{
	TITLE("Average","Routine");

	cal();

	return;
}


void Average::cal()
{
	TITLE("Average","cal");
	cout << " Calculate the averaged ion positions." << endl;
	cout << " Geometry ranges from " << INPUT.geo_1 << " to " << INPUT.geo_2 << endl;
	cout << " atom number is " << INPUT.natom << endl;
	cout << " geo_interval is " << INPUT.geo_interval << endl;
	
	double cella, cellb, cellc;
	double* xatom = new double[INPUT.natom];
	double* yatom = new double[INPUT.natom];
	double* zatom = new double[INPUT.natom];
	double* xvel = new double[INPUT.natom];
	double* yvel = new double[INPUT.natom];
	double* zvel = new double[INPUT.natom];

	string* element = new string[INPUT.natom];
	for(int i=0; i<INPUT.natom; ++i)
	{
		xatom[i] = 0.0;
		yatom[i] = 0.0;
		zatom[i] = 0.0;
		xvel[i] = 0.0;
		yvel[i] = 0.0;
		zvel[i] = 0.0;
	}

	int count_geometry_number = 0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
	{
		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) continue;

		CellFile cel;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
		++count_geometry_number;

		// calculate the ionic density
		const double rho_ion = INPUT.natom / cel.volume;
		if(count_geometry_number==1)
		{
			cella = cel.a1.x;
			cellb = cel.a2.y;
			cellc = cel.a3.z;

			cout << " Volume of the input cell = " << cel.volume << " A^3" << endl;
			cout << " Average ion density = " << rho_ion << endl;
			// Fermi vector
			double kf = pow(3*PI*PI*rho_ion,1.0/3.0);
	        cout << " Fermi vector = " << kf << endl; 
	        assert(kf > 0.0);
		}


		// read in average, we constrcut a new CellFile
		// because when reading the velocity,
		// all things in CellFile are generated again.
		// If we use the above one, segmental fault will happen
		CellFile cel2;
		cel2.file_name = ss.str();

		if( !CellFile::ReadVelocity( cel2 ) )
		{
			cout << " Strange! Can't find the velocity file!" << endl;
			exit(0);
		}



		for(int i=0; i<INPUT.natom; ++i)
		{
			element[i] = cel.atom[0].id;
			xatom[i] += cel.atom[0].pos[i].x;
			yatom[i] += cel.atom[0].pos[i].y;
			zatom[i] += cel.atom[0].pos[i].z;
			xvel[i] += cel2.atom[0].vel[i].x;
			yvel[i] += cel2.atom[0].vel[i].y;
			zvel[i] += cel2.atom[0].vel[i].z;
		}
	}

	// get the averaged atom positions
	for(int i=0; i<INPUT.natom; ++i)
	{
		xatom[i] = xatom[i]/count_geometry_number;
		yatom[i] = yatom[i]/count_geometry_number;
		zatom[i] = zatom[i]/count_geometry_number;
		xvel[i] = xvel[i]/count_geometry_number;
		yvel[i] = yvel[i]/count_geometry_number;
		zvel[i] = zvel[i]/count_geometry_number;
	}


	cout << " movement(Angstrom)  input by user: ";
	cout << INPUT.movement_x << " " << INPUT.movement_y << " " << INPUT.movement_z << endl;
	cout << " this movement is used to avoid creating vacancy." << endl;
	for(int i=0; i<INPUT.natom; ++i)
	{
		xatom[i] += INPUT.movement_x;
		yatom[i] += INPUT.movement_y;
		zatom[i] += INPUT.movement_z;
	}


	// confine all the atoms within the cell
	for(int i=0; i<INPUT.natom; ++i)
	{
		while(xatom[i]<0) xatom[i] += cella;
		while(yatom[i]<0) yatom[i] += cellb;
		while(zatom[i]<0) zatom[i] += cellc;

		while(xatom[i]>=cella) xatom[i] -= cella;
		while(yatom[i]>=cellb) yatom[i] -= cellb;
		while(zatom[i]>=cellc) zatom[i] -= cellc;
	}


	// output the final pair distribution function
	ofstream ofs(INPUT.geo_out.c_str());
	cout << " print out the final geometry into file: " << INPUT.geo_out << endl;
	ofs << "%BLOCK LATTICE_CART" << endl;
	ofs << cella << " 0 0" << endl;
	ofs << "0 " << cellb << " 0" << endl;
	ofs << "0 0 " << cellc << endl;
	ofs << "%ENDBLOCK LATTICE_CART" << endl;
	ofs << "%BLOCK POSITIONS_CART" << endl;

	// print out the geometry into .xyz format
	ofstream xyz("plot.xyz");
	cout << " print out the geometry into .xyz file: " << "plot.xyz" << endl;
	xyz << INPUT.natom << endl;
	xyz << "ForPloting" << endl;

	for(int i=0; i<INPUT.natom; ++i)
	{
		ofs << "Li" << " " << xatom[i] << " " << yatom[i] << " " << zatom[i] << endl;
		xyz << "Li" << " " << xatom[i] << " " << yatom[i] << " " << zatom[i] << endl;
	}
	
	ofs << "%END BLOCK POSITIONS_CART" << endl;
	ofs.close();
	xyz.close();

	// output the velocity
	cout << " print out the final velocity info into file: " << INPUT.vel_out << endl;
	ofstream ofsv(INPUT.vel_out.c_str());
	ofsv << setprecision(12);
	ofsv << " MD STEP= 1" << endl;
	ofsv << " 0 <===xLogS (a.u.)" << endl;
	ofsv << " 0 <===vLogS (a.u.)" << endl;
	ofsv << " 0 <=== barostat velocity 1" << endl;
	ofsv << " 0 <=== barostat velocity 2" << endl;
	ofsv << " 0 <=== barostat velocity 3" << endl;
	ofsv << "ION VELOCITIES (a.u.):" << endl;
	for(int i=0; i<INPUT.natom; ++i)
	{
		ofsv << xvel[i] << " " << yvel[i] << " " << zvel[i] << endl;
	}
	ofsv.close();


	// delete stuff
	delete[] xatom;
	delete[] yatom;
	delete[] zatom;
	delete[] element; 

	return;
}
