#include "cellFile.h"
#include "input.h"
#include "mdp3.h"
#include "mdp.h"
#include "math.h"
#include "water.h"
#include "HBs.h"

MDP3::MDP3()
{
}

MDP3::~MDP3(){}

void MDP3::Routine()
{
	TITLE("MDP3","Routine");
	
	ofs_running << "Analysis based on Mean Density Profile and Instantaneous Liquid Interface (ILI)" << endl;
	ofs_running << "Recommend paper: J. Phys. Chem. B 2010, 114, 1954-1958." << endl;

    this->nx = INPUT.u1;
    this->ny = INPUT.u2;
	this->x0 = INPUT.x0; 
	this->y0 = INPUT.y0; 
	this->dx = INPUT.dx; 
	this->dy = INPUT.dy; 

	this->coord_xy = new double*[nx];
    for(int ix=0; ix<nx; ++ix)
    {
        this->coord_xy[ix] = new double[ny]();
    }
	
	// setup interface for each snapshot and its
	// associated gradient
	this->interface = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->interface[ix] = new double[INPUT.ny]();
	}
	this->gradient = new double**[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->gradient[ix] = new double*[INPUT.ny];
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			this->gradient[ix][iy] = new double[3]();
		}
	}

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	// input ili file
	ifstream ifs(INPUT.ili_file.c_str());
	if(!ifs)
	{
		cout << "Cannot find the ILI file: " << INPUT.ili_file << endl;
		exit(0);
	}
	else
	{
		cout << "Open the ILI file: " << INPUT.ili_file << endl;
	}
	

//-----------------  CORE CODE -------------------------------------

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "snapshot" << setw(12) << igeo << endl;

		proximity(ifs, cel);
	}	

//------------------------------------------------------------------
	// print out the information
	ofstream ofs2D("dis_2D.dat");

    double sum=0.0;
    for(int iy=0; iy<ny; ++iy)
    {
        for(int ix=0; ix<nx; ++ix)
        {
            sum += coord_xy[ix][iy] * INPUT.dx * INPUT.dy;
        }
    }

	cout << "sum = " << sum << endl;
	ofs_running << "sum = " << sum << endl;

    if(sum>0.0)
    {
        for(int iy=0; iy<ny; ++iy)
        {
            for(int ix=0; ix<nx; ++ix)
            {
                ofs2D << coord_xy[ix][iy]/sum << " ";
            }
            ofs2D << endl;
        }
    }
	
	//-------------
	// clean up
	//-------------
	// close files
	ifs.close();
	ofs2D.close();

	// delete arrays
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		delete[] interface[ix];
	}
	delete[] interface;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			delete[] gradient[ix][iy];
		}
		delete[] gradient[ix];
	}
	delete[] gradient;

    for(int ix=0; ix<nx; ++ix)
    {
        delete[] coord_xy[ix];
    }
    delete[] coord_xy;

	cout << "All done. Have a great day!" << endl;
	ofs_running << "All done. Have a great day!" << endl;

	return;
}


void MDP3::proximity(ifstream &ifs, const Cell &cel)
{
	TITLE("MDP","proximity");

	// get ito, ith, and itc.
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

	// read the Instantaneous Liquid Interface (ILI)
	MDP::read_ili(ifs, this->interface, this->gradient);

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);

	// print out the distance of each atom to the ILI
	for(int ia=0; ia<cel.atom[ito].na;++ia)
	{
		if(water[ia].nH!=2) continue;
		if(water[ia].ndonate==2) continue; // mohan added 2018-08-25


		double posx = cel.atom[ito].pos[ia].x;
		double posy = cel.atom[ito].pos[ia].y;
		double posz = cel.atom[ito].pos[ia].z;
		while( posx >= norm1 ) posx -= norm1;	
		while( posx < 0 ) posx += norm1;
		while( posy >= norm2 ) posy -= norm2;	
		while( posy < 0 ) posy += norm2;
		if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide") 
		{
			while( posz >= norm3 ) posz -= norm3;	
			while( posz <0 ) posz += norm3;
		}
		else while( posz > INPUT.upper_z ) posz -= norm3;	


		int six=-1; int siy=-1;
		double drx=0.0; double dry=0.0; double drz=0.0;
		MDP::which_surface(six, siy, norm1, norm2, posx, posy, posz, drx, dry, drz, this->interface);
		double aaa = -(drx * gradient[six][siy][0] + dry * gradient[six][siy][1] + drz * gradient[six][siy][2]);
		//cout << "posz=" << posz << " aaa=" << aaa << endl;

		int indexX = (aaa-x0)/dx;

		for(int ii=0; ii<2; ++ii)
		{
			int ih=water[ia].indexH[ii];
			double ho_dx = shortest(posx, cel.atom[ith].pos[ih].x, INPUT.celldm1);
			double ho_dy = shortest(posy, cel.atom[ith].pos[ih].y, INPUT.celldm2);
			double ho_dz = shortest(posz, cel.atom[ith].pos[ih].z, INPUT.celldm3);
			double ho_dis = sqrt(ho_dx*ho_dx + ho_dy*ho_dy + ho_dz*ho_dz);
			double ucos1=(ho_dx*gradient[six][siy][0]+ho_dy*gradient[six][siy][1]+ho_dz*gradient[six][siy][2])/ho_dis;
			double degree=acos(ucos1)/3.1415926*180;
			//cout << aaa << " " << acos(ucos1)/3.1415926*180 << endl;
			//cout << acos(1) << endl;

			int indexY = (degree-y0)/dy;

			if(indexX<nx and indexY<ny and indexX>=0 and indexY>=0)
			{
				coord_xy[indexX][indexY]+=1.0;
			}
		}

	}// end ia

	delete[] water;

	return;
}
