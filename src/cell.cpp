#include "cell.h"
#include "matrix3.h"
Cell::Cell()
{
	read_and_used=false;
	nat=0;	
	nbonds=0;
	nangles=0;
	ntype=0;
	snapshot_index = -1;
	snapshot_time = -1.0;
}

Cell::~Cell()
{
	delete[] atom;
	delete[] wan_centers;
	delete[] eig;
}

double Cell::last_volume = 0.0;

// transform the coordinates from direct to cartesian
void Cell::direct2cartesian(const int &it, const int &i)
{
	assert( atom[it].pos != nullptr);
	assert( atom[it].posd != nullptr);
	atom[it].pos[i].x = atom[it].posd[i].x * a1.x + atom[it].posd[i].y * a2.x + atom[it].posd[i].z * a3.x;
	atom[it].pos[i].y = atom[it].posd[i].x * a1.y + atom[it].posd[i].y * a2.y + atom[it].posd[i].z * a3.y;
	atom[it].pos[i].z = atom[it].posd[i].x * a1.z + atom[it].posd[i].y * a2.z + atom[it].posd[i].z * a3.z;
	//cout << " cartesian : " << atom[it].pos[i].x << " " << atom[it].pos[i].y << " " << atom[it].pos[i].z << endl;
	return;
}	

// transform the coordinates from cartesian to direct 
void Cell::cartesian2direct(const int &it, const int &i)
{
	assert( atom[it].pos != nullptr);
	assert( atom[it].posd != nullptr);
	
	Matrix3 lattice_vector, inv_lat;
	lattice_vector.e11=a1.x;
	lattice_vector.e12=a1.y;
	lattice_vector.e13=a1.z;
	lattice_vector.e21=a2.x;
	lattice_vector.e22=a2.y;
	lattice_vector.e23=a2.z;
	lattice_vector.e31=a3.x;
	lattice_vector.e32=a3.y;
	lattice_vector.e33=a3.z;

	inv_lat = lattice_vector.Inverse();
	Vector3<double> direct_vec, cartesian_vec;
	cartesian_vec.x = atom[it].pos[i].x;
	cartesian_vec.y = atom[it].pos[i].y;
	cartesian_vec.z = atom[it].pos[i].z;
	direct_vec = cartesian_vec * inv_lat;
	atom[it].posd[i].x = direct_vec.x;
	atom[it].posd[i].y = direct_vec.y;
	atom[it].posd[i].z = direct_vec.z;
	return;
}	



void Cell::add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		double &x, double &y, double &z) const
{
	x = this->atom[it].pos[ia].x + i*a1.x + j*a2.x + k*a3.x;
	y = this->atom[it].pos[ia].y + i*a1.y + j*a2.y + k*a3.y;
	z = this->atom[it].pos[ia].z + i*a1.z + j*a2.z + k*a3.z;
	return;
}

void Cell::add_cell_length(const int &it, const int &ia, const int &i, const int &j, const int &k, 
		float &x, float &y, float &z) const
{
	x = this->atom[it].pos[ia].x + i*a1.x + j*a2.x + k*a3.x;
	y = this->atom[it].pos[ia].y + i*a1.y + j*a2.y + k*a3.y;
	z = this->atom[it].pos[ia].z + i*a1.z + j*a2.z + k*a3.z;
	return;
}


void Cell::Print_water(ofstream &ofs, const int &ito, const int &ith, 
const double &norm1, const double &norm2, const double &norm3, const double &rcut_oh,
const int &output_type)
{
	Water *water = new Water[atom[ito].na];

	ofs << endl;
	if(output_type==1)
	{
		ofs << "Atoms" << endl;
		ofs << endl;
	}
	else if(output_type==2)
	{
		ofs << "Bonds" << endl;
		ofs << endl;
	}
	else if(output_type==3)
	{
		ofs << "Angles" << endl;
		ofs << endl;
	}

	// 1) the first part: search for O-H pair
	int iat=0;
	int cc=0;
	int iwater=0;
	int ibond=0;
	for(int ia=0; ia<atom[ito].na; ++ia)
	{
		for(int ia2=0; ia2<atom[ith].na; ++ia2)
		{
			double dis = distance(atom[ito].pos[ia], atom[ith].pos[ia2], norm1, norm2, norm3);
			if(dis < rcut_oh)
			{
				int ind=water[ia].nH;
				water[ia].indexH[ind]=ia2;
				water[ia].disH[ind]=dis;
				water[ia].nH++;	
			}
		}
		if(water[ia].nH!=2) cout << ia+1 << " O atom, #HBs: " << water[ia].nH << endl;
		else  ++cc;


		if(output_type==1) // water coordinates
		{
			double ox = atom[ito].pos[ia].x;
			double oy = atom[ito].pos[ia].y;
			double oz = atom[ito].pos[ia].z;
			++iat;
			++iwater;	
			ofs << iat << " "
				<< iwater << " " 
				<< ito+1 << " " 
				<< atom[ito].charge << " "  
				<< ox << " " << oy << " " << oz << endl;
			for(int ia2=0; ia2<water[ia].nH; ++ia2)
			{
				const int ind = water[ia].indexH[ia2];
				double dx = shortest(ox, atom[ith].pos[ind].x, norm1); 
				double dy = shortest(oy, atom[ith].pos[ind].y, norm2); 
				double dz = shortest(oz, atom[ith].pos[ind].z, norm3); 
				++iat;
				ofs << iat << " " 
					<< iwater << " " 
					<< ith+1 << " " 
					<< atom[ith].charge << " " 
					<< ox-dx << " " 
					<< oy-dy << " " 
					<< oz-dz << endl;
			}
		}
		else if(output_type==2) // water bonds
		{
			++iat;
			int iao = iat;
			for(int ia2=0; ia2<water[ia].nH; ++ia2)
			{ 
				++ibond;
				++iat;
				const int ind = water[ia].indexH[ia2];
				ofs << ibond << " 1 " 
					<< iao << " " << iat << endl;
			}
		}
		else if(output_type==3) // water angles between atoms
		{
			ofs << ia+1 << " 1 " << iat+2 << " " << iat+1 << " " << iat+3 << endl;	
			++iat;
			++iat;
			++iat;
		}
	}

	cout << "2 HBs, no problem " << cc << endl;
	delete[] water;
	return;
}


void Cell::read_wannier_centers(ifstream &ifs_wan, const int &nbands)
{
//	cout << "reading wannier centers now." << endl;

	int frame; double time;
	string useless;
	ifs_wan >> frame;
	cout << "snapshot " << frame << " (wannier)" << endl;
	READ_VALUE(ifs_wan, useless);
	time = frame*INPUT.dt_snapshots;
	assert(this->snapshot_index == frame);
//	ofs_running << frame << " " << time << endl;

	double wx, wy, wz;
	for(int i=0; i<nbands; ++i)
	{
		ifs_wan >> wx >> wy >> wz;
		this->wan_centers[i].x = wx*0.529177; 
		this->wan_centers[i].y = wy*0.529177;
		this->wan_centers[i].z = wz*0.529177; // unit from bohr to Angstroms
//		ofs_running << "read " << i << " " << wx << " " << wy << " " << wz << endl; 
	}
	return;
}


// mohan added 2017-04-08
void Cell::read_pos_ili(ifstream &ifs_pos_ili, const int &it)
{
	//cout << "read_pos_ili " << atom[it].na << endl; 
	string id;
	double x,y,z;
	double dis=0.0;
	for(int ia=0; ia< this->atom[it].na; ++ia)
	{
		ifs_pos_ili >> id >> x >> y >> z >> dis;
		atom[it].pos_ili[ia] = dis;
		//cout << it+1 << " " << ia+1 << " distance to ILI (Angstroms) is " <<  dis << endl;
	}

	return;
}


void Cell::read_eig(ifstream &ifs_eig, const int &nbands)
{
//	cout << "Read eigenvalues." << endl;
	string line1;
	string line2;
	READ_VALUE(ifs_eig, line1);
	READ_VALUE(ifs_eig, line2);	
//	cout << line1 << endl;
	
	for(int i=0; i<nbands; ++i)
	{
		ifs_eig >> this->eig[i];
	}

	return;
}

void Cell::clean()
{
	delete[] this->atom;
	delete[] this->wan_centers;
	delete[] this->eig;
	this->atom = nullptr;
	this->wan_centers = nullptr;
	this->eig = nullptr;
	return;
}

void Cell::atom_mass()
{
	if(this->atom == nullptr) return;
	for(int it=0; it<this->ntype; ++it)
	{
		this->atom[it].cal_mass();
	}
}

void Cell::init_cel(Input& Inp)
{
	this->ntype = Inp.ntype;
	this->nat = Inp.natom;
	delete[] this->atom; 
	this->atom = new Atoms[ntype];
	if(ntype==1)
	{
		this->atom[0].na = this->nat;
		this->atom[0].id = INPUT.id1;
	}
	else 
	{
		this->atom[0].na = INPUT.natom1; this->atom[0].id = INPUT.id1;
		if(INPUT.ntype>=2) {this->atom[1].na = INPUT.natom2; this->atom[1].id = INPUT.id2; }
		if(INPUT.ntype>=3) {this->atom[2].na = INPUT.natom3; this->atom[2].id = INPUT.id3; }
		if(INPUT.ntype>=4) {this->atom[3].na = INPUT.natom4; this->atom[3].id = INPUT.id4; }
		if(INPUT.ntype>=5) 
		{
			cout<<"We do not support ntype > 4 yet."<<endl;
			exit(0);
		}
	}
	for(int it=0; it<ntype; ++it)
    {
        delete[] this->atom[it].pos;	this->atom[it].pos = new Vector3<double>[this->atom[it].na];
		delete[] this->atom[it].posd;   this->atom[it].posd = new Vector3<double>[this->atom[it].na];
		if(INPUT.read_velocity)
		{
			delete[] this->atom[it].vel; this->atom[it].vel = new Vector3<double>[this->atom[it].na];
		}
    }
	return;
}

double Cell::cal_volume()
{
	this->volume = this->a1.x*this->a2.y*this->a3.z + this->a1.y*this->a2.z*this->a3.x + this->a1.z*this->a2.x*this->a3.y -
	  this->a1.x*this->a2.z*this->a3.y - this->a1.y*this->a2.x*this->a3.z - this->a1.z*this->a2.y*this->a3.x;
	this->last_volume = this->volume;
	return this->volume;
}