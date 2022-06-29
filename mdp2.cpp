#include "cellFile.h"
#include "input.h"
#include "mdp.h"
#include "mdp2.h"
#include "math.h"
#include "water.h"
#include "HBs.h"

void MDP2::Routine()
{
	TITLE("MDP2","Routine");
	
	cout << " === Compute Mean Density Profile based on Mean Liquid Interface ===" << endl;

	// (1) set up the basic parameters: delta r in real space.
	this->dr = INPUT.dr;

	// radius cutoff in real space, usually I choose a/2, where a is the lattice constant.
	this->rcut = INPUT.rcut;

	// number of radial mesh.
	this->nmesh = int(rcut / dr) +  1;

	this->mean_density = new double[nmesh]();
	for(int i=0; i<nmesh; ++i)
	{
		mean_density[i]=0.0;
	}

	// (2) setup interface and its gradient arrays
	this->interface = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->interface[ix] = new double[INPUT.ny];
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			this->interface[ix][iy] = 0.0;
		}
	}
	this->gradient = new double**[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->gradient[ix] = new double*[INPUT.ny];
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			this->gradient[ix][iy] = new double[3];
			this->gradient[ix][iy][0] = 0.0; 
			this->gradient[ix][iy][1] = 0.0; 
			this->gradient[ix][iy][2] = 0.0; 
		}
	}


	// (3) setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	//cout << "geo1=" << INPUT.geo_1 << endl;
	//cout << "geo2=" << INPUT.geo_2 << endl;

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

	// read in the averaged surface first.
	if(INPUT.func==2)
	{
		read_write_ili(ifs);
	}

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
		cout << "igeo=" << igeo << endl;

		// last parameter 1: read and write
		// 2: compute average mean density profile
		proximity(ifs, cel, INPUT.func);

	}	

	// output data
	string filename;
	if(INPUT.func==1)
	{
		filename = "ave_ili.dat";
	}
	else if(INPUT.func==2)
	{
		filename = "ave_mdp_result.dat";
	}
	else
	{
		cout << "warning for INPUT.func." << endl;
		exit(0);
	}
	ofstream ofs(filename.c_str());

	if(INPUT.func==2)
	{
		ofs << "x O Oacc Odon H" << endl;

		double ave_na = 0.0;
		for(int ir=0; ir<nmesh; ++ir)
		{
			ave_na += mean_density[ir];
		}
		ave_na/=(double)count_geometry_number;
		cout << "ave_na = " << ave_na << endl;
		cout << "celldm1= " << INPUT.celldm1 << endl;
		cout << "celldm2= " << INPUT.celldm2 << endl;
		cout << "dr= " << dr << endl;

		// 96 water molecules in a 12.445*12.445*18.5162 cell 
		// has a density of 1.0 g/cm3
		// 96/12.445/12.445/18.5162 = .0334756923 number_of_O/Angstrom^3

		const double fac = (double)count_geometry_number*(INPUT.celldm1*INPUT.celldm2*dr)*0.03347569;
		assert(fac!=0.0);
		cout << fac << endl;
		for(int ir=0; ir<nmesh; ++ir)
		{
			mean_density[ir] /= fac; 
			ofs << ir*dr-2.0 << " " << mean_density[ir] << endl; 
		}
	}
	else if(INPUT.func==1)
	{	
		if(RANK==0)
		{
			assert(count_geometry_number>0);
			ofs << "ILI_AVERAGED" << endl;
			ofs << INPUT.nx << " " << INPUT.ny << endl;
			for(int ix=0; ix<INPUT.nx; ++ix)
			{
				for(int iy=0; iy<INPUT.ny; ++iy)
				{
					interface[ix][iy]/=(double)count_geometry_number;
					ofs << interface[ix][iy] << " ";
				}
				ofs << endl;
			}
			for(int j=0; j<3; ++j)
			{
				ofs << "GRAD" << j+1 << endl;
				ofs << INPUT.nx << " " << INPUT.ny << endl;
				for(int ix=0; ix<INPUT.nx; ++ix)
				{
					for(int iy=0; iy<INPUT.ny; ++iy)
					{
						gradient[ix][iy][j]/=(double)count_geometry_number;
						ofs << gradient[ix][iy][j] << " ";
					}
					ofs << endl;
				}
			}
		}
	}

	// close files
	ifs.close();
	ofs.close();

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

	// delete
	delete[] mean_density;
	return;
}


void MDP2::proximity(ifstream &ifs, const Cell &cel, const int &func)
{
	TITLE("MDP2","proximity");

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

	if(func==1)
	{
		read_write_ili(ifs);
		return;
	}


	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];

	// search for OH
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
		{
			double ho_dis = distance(cel.atom[ito].pos[ia],cel.atom[ith].pos[ia2],norm1,norm2,norm3);
			if(ho_dis<INPUT.rcut_oh)
			{
				int ind=water[ia].nH;
				water[ia].indexH[ind]=ia2;
				water[ia].disH[ind]=ho_dis;
				water[ia].nH++;	
			}
		}
	}


	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		double posx = cel.atom[ito].pos[ia].x;
		double posy = cel.atom[ito].pos[ia].y;
		double posz = cel.atom[ito].pos[ia].z;
		while( posx >= norm1 ) posx -= norm1;	
		while( posx < 0 ) posx += norm1;
		while( posy >= norm2 ) posy -= norm2;	
		while( posy < 0 ) posy += norm2;
		if(INPUT.system=="water")
		{
			while( posz >= norm3 ) posz -= norm3;	
		}
		else
		{
			while( posz > INPUT.upper_z ) posz -= norm3;	
		}

		int six=-1;
		int siy=-1;
		double drx=0.0;
		double dry=0.0;
		double drz=0.0;
		MDP::which_surface(six, siy, norm1, norm2, posx, posy, posz, drx, dry, drz, this->interface);

		double aaa = drx * gradient[six][siy][0] + dry * gradient[six][siy][1] + drz * gradient[six][siy][2];
		aaa=-aaa;

		if(aaa>=-2)
		{
			//int index = int(aaa/this->dr);
			int index = int((aaa+2.0)/this->dr);
			if(index>=0 and index<this->nmesh)
			{
				mean_density[index]+=1.0;
			}
		}
		else
		{
			cout << setw(10) << ia << setw(15) << aaa 
				<< " pos " << posx << " " << posy << " " << posz 
				<< " dr " << drx << " " << dry << " " << drz 
				<< " grad " << gradient[six][siy][0] << " " << gradient[six][siy][1] << " " << gradient[six][siy][2] << endl;
		}
	}

	delete[] water;

	return;
}

void MDP2::read_write_ili(ifstream &ifs)
{
	string title;
	READ_VALUE(ifs, title);
	//cout << title << endl;
	int nx, ny;
	ifs >> nx >> ny;
	assert(INPUT.nx == nx);
	assert(INPUT.ny == ny);

	double tmp=0.0;
	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			ifs >> tmp; 
			interface[ix][iy] += tmp;
		}
	}

	double*** tmpg = new double**[nx]; 
	for(int ix=0; ix<nx; ++ix)
	{
		tmpg[ix] = new double*[ny];
		for(int iy=0; iy<ny; ++iy)
		{
			tmpg[ix][iy] = new double[3];
			tmpg[ix][iy][0] = 0.0; 
			tmpg[ix][iy][1] = 0.0; 
			tmpg[ix][iy][2] = 0.0; 
		}
	}

	for(int j=0; j<3; ++j)
	{
		READ_VALUE(ifs,title);
		//cout << title << endl;
		ifs >> nx >> ny;
		for(int ix=0; ix<nx; ++ix)
		{
			for(int iy=0; iy<ny; ++iy)
			{
				ifs >> tmpg[ix][iy][j]; 
			}
		}
	}

	// unit 
	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			double sum = sqrt(tmpg[ix][iy][0]*tmpg[ix][iy][0]+
			tmpg[ix][iy][1]*tmpg[ix][iy][1]+
			tmpg[ix][iy][2]*tmpg[ix][iy][2]);
	
			assert(sum>0.0);
				
			for(int j=0; j<3; ++j)
			{
				tmpg[ix][iy][j]/=sum;
				gradient[ix][iy][j]+=tmpg[ix][iy][j];
			}
		}
	}


	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			delete[] tmpg[ix][iy];
		}
		delete[] tmpg[ix];
	}
	delete[] tmpg;



	return;
}
