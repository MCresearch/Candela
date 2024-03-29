#include "atoms.h"
#include "input.h"
#include "const.h"

Atoms::Atoms()
{
	na=0;
	mass=0;
	pos = nullptr;
	posd = nullptr;
	vel = nullptr;
	pos_ili = nullptr;
}


Atoms::~Atoms()
{
	//cout << "deallocate atoms" << endl;
	delete[] pos;
	delete[] posd;
	delete[] vel;
	delete[] pos_ili;
}

// read the 'fractional' coordinates of atoms
// or the 'cartesian' coordinates of atoms.
void Atoms::read_pos(ifstream &ifs,bool frac)
{
	assert(na>0);
	string idtmp;
	for(int i=0; i<na; ++i)
	{
		ifs >> idtmp;
		if(frac) // fractional coordinates
		{
			ifs >> posd[i].x >> posd[i].y >> posd[i].z;
			// mohan add 2016-03-05
			while(posd[i].x<0){posd[i].x+=1.0;}
			while(posd[i].x>1){posd[i].x-=1.0;}
			while(posd[i].y<0){posd[i].y+=1.0;}
			while(posd[i].y>1){posd[i].y-=1.0;}
			while(posd[i].z<0){posd[i].z+=1.0;}
			while(posd[i].z>1){posd[i].z-=1.0;}
			//cout << " pos=" << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl; 
		//	cout << ifs.good() << " " << ifs.bad() << " " << ifs.fail() << endl;
		}	
		else // cartesian coordinates
		{
			ifs >> pos[i].x >> pos[i].y >> pos[i].z;

			// right now only works for cartesian coordinates 2015-06-22
		    if(INPUT.movement_x != 0 || INPUT.movement_y !=0 || INPUT.movement_z !=0)
			{  
				pos[i].x += INPUT.movement_x;	
				pos[i].y += INPUT.movement_y;	
				pos[i].z += INPUT.movement_z;	
			}
			//cout << " pos=" << pos[i].x << " " << pos[i].y << " " << pos[i].z << endl;
		}
		if( ifs.fail() )
		{
			cout << " Reading atom " << i+1 << endl;
			QUIT("Error in reading atom positions");
		}
	}
	return;
}


// read the 'fractional' coordinates of atoms
// or the 'cartesian' coordinates of atoms.
void Atoms::read_pos_2(ifstream &ifs,bool frac)
{
	assert(na>0);
	string tmp;
	for(int i=0; i<na; ++i)
	{
		if(frac)
		{
			ifs >> posd[i].x >> posd[i].y >> posd[i].z;
			getline(ifs, tmp);
			//cout << " posd=" << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl; 

		//	cout << ifs.good() << " " << ifs.bad() << " " << ifs.fail() << endl;
			if( ifs.fail() )
			{
				cout << " Reading atom " << i+1 << endl;
				QUIT("Error in reading atom positions 2");
			}
		}	
		else
		{
			ifs >> pos[i].x >> pos[i].y >> pos[i].z;
			getline(ifs, tmp);
		}
	}
	return;
}



void Atoms::read_pos_3(ifstream &ifs)
{
    assert(na>0);
    string tmp;
    for(int i=0; i<na; ++i)
    {
        ifs >> id >> posd[i].x >> posd[i].y >> posd[i].z;
        getline(ifs, tmp);
//		cout << id << " " << posd[i].x << " " << posd[i].y << " " << posd[i].z << endl;
    }
    return;
}


void Atoms::read_pos_4(ifstream &ifs, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &a3)
{
	assert(na>0);
	string tmp;
	for(int i=0; i<na; ++i)
	{
		ifs >> pos[i].x >> pos[i].y;
		READ_VALUE(ifs, pos[i].z);
		//cout << pos[i].x << " " << pos[i].y << " " << pos[i].z << endl;
		// from bohr to Angstrom
		pos[i].x *= BOHR;
		pos[i].y *= BOHR;
		pos[i].z *= BOHR;

		// mohan updated 2018-04-26
		// correct
		while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a1.y; pos[i].z -= a1.z; }
		while( pos[i].x < 0){ pos[i].x += a1.x; pos[i].y += a1.y; pos[i].z += a1.z; }
		while( pos[i].y >= a2.y ){ pos[i].x -= a2.x; pos[i].y -= a2.y; pos[i].z -= a2.z; }
		while( pos[i].y < 0){ pos[i].x += a2.x; pos[i].y += a2.y; pos[i].z += a2.z; }
		while( pos[i].z >= a3.z ){ pos[i].x -= a3.x; pos[i].y -= a3.y; pos[i].z -= a3.z; }
		while( pos[i].z < 0){ pos[i].x += a3.x; pos[i].y += a3.y; pos[i].z += a3.z; }

	}
	return;
}

void Atoms::read_pos_5(ifstream &ifs, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &a3, double &lat_const)
{
	assert(na>0);
	string tmp;
	for(int i=0; i<na; ++i)
	{
		ifs >> tmp >> tmp;

		ifs >> pos[i].x >> pos[i].y;
		READ_VALUE(ifs, pos[i].z);
		//cout << pos[i].x << " " << pos[i].y << " " << pos[i].z << endl;
		// from bohr to Angstrom
		pos[i].x *= BOHR*lat_const;
		pos[i].y *= BOHR*lat_const;
		pos[i].z *= BOHR*lat_const;

		// mohan updated 2018-04-26
		// correct
		while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a1.y; pos[i].z -= a1.z; }
		while( pos[i].x < 0){ pos[i].x += a1.x; pos[i].y += a1.y; pos[i].z += a1.z; }
		while( pos[i].y >= a2.y ){ pos[i].x -= a2.x; pos[i].y -= a2.y; pos[i].z -= a2.z; }
		while( pos[i].y < 0){ pos[i].x += a2.x; pos[i].y += a2.y; pos[i].z += a2.z; }
		while( pos[i].z >= a3.z ){ pos[i].x -= a3.x; pos[i].y -= a3.y; pos[i].z -= a3.z; }
		while( pos[i].z < 0){ pos[i].x += a3.x; pos[i].y += a3.y; pos[i].z += a3.z; }

	}
	return;
}

void Atoms::read_vel(ifstream &ifs)
{
	assert(na>0);
//	cout << " atom number is " << na << endl;
	for(int i=0; i<na; ++i)
	{
		ifs >> vel[i].x >> vel[i].y >> vel[i].z;
		//cout << " vel=" << vel[i].x << " " << vel[i].y << " " << vel[i].z << endl;
	}
	return;
}

void Atoms::cal_mass()
{
	if(this->id == "H")
	{
		this->mass = MASS_H;
	}
	else if(this->id == "D")
	{
		this->mass = MASS_D;
	}
	else if(this->id == "T")
	{
		this->mass = MASS_T;
	}
	else if(this->id == "Li")
	{
		this->mass = MASS_Li;
	}
	else if(this->id == "Be")
	{
		this->mass = MASS_Be;
	}
	else if(this->id == "B")
	{
		this->mass = MASS_B;
	}
	else if(this->id == "C")
	{
		this->mass = MASS_C;
	}
	else if(this->id == "O")
	{
		this->mass = MASS_O;
	}
	else if(this->id == "F")
	{
		this->mass = MASS_F;
	}
	else if(this->id == "Na")
	{
		this->mass = MASS_Na;
	}
	else if(this->id == "Mg")
	{
		this->mass = MASS_Mg;
	}
	else if(this->id == "Al")
	{
		this->mass = MASS_Al;
	}
	else if(this->id == "Si")
	{
		this->mass = MASS_Si;
	}
	else
	{
		// cout<<"Warning: we do not support the mass of "<<this->id<<" yet!"<<endl;
	}
}


