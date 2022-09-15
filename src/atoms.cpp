#include "atoms.h"
#include "input.h"

Atoms::Atoms()
{
	na=0;
	mass=0;
	allocate_pos=false;
	allocate_posd=false;
	allocate_vel=false;
	allocate_pos_ili=false;
	pos = nullptr;
	posd = nullptr;
	vel = nullptr;
}


Atoms::~Atoms()
{
	//cout << "deallocate atoms" << endl;
	if(allocate_pos==true) delete[] pos;
	if(allocate_posd==true) delete[] posd;
	if(allocate_vel==true) delete[] vel;
	if(allocate_pos_ili==true) delete[] pos_ili;
}

// read the 'fractional' coordinates of atoms
// or the 'cartesian' coordinates of atoms.
void Atoms::read_pos(ifstream &ifs,bool frac)
{
	assert(na>0);
	string idtmp;
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	allocate_pos = true;
	allocate_posd = true;
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
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	allocate_pos = true;
	allocate_posd = true;
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
    pos = new Vector3<double>[na];
    posd = new Vector3<double>[na];
	allocate_pos = true;
	allocate_posd = true;
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
	const double norm1 = INPUT.celldm1; 
	const double norm2 = INPUT.celldm2; 
	const double norm3 = INPUT.celldm3;
	assert(na>0);
	string tmp;
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	allocate_pos = true;
	allocate_posd = true;
	for(int i=0; i<na; ++i)
	{
		ifs >> pos[i].x >> pos[i].y;
		READ_VALUE(ifs, pos[i].z);
		//cout << pos[i].x << " " << pos[i].y << " " << pos[i].z << endl;
		// from bohr to Angstrom
		pos[i].x *= BOHR;
		pos[i].y *= BOHR;
		pos[i].z *= BOHR;

		// mohan added 2016-11-09
	//	while( pos[i].x >= norm1 ) pos[i].x -= norm1;	
	//	while( pos[i].x < 0 ) pos[i].x += norm1;
	//	while( pos[i].y >= norm2 ) pos[i].y -= norm2;	
	//	while( pos[i].y < 0 ) pos[i].y += norm2;
	//	while( pos[i].z >= norm3 ) pos[i].z -= norm3;	
	//	while( pos[i].z < 0 ) pos[i].z += norm3;
	
		// mohan updated 2018-04-26
		// correct
		while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a1.y; pos[i].z -= a1.z; }
		while( pos[i].x < 0){ pos[i].x += a1.x; pos[i].y += a1.y; pos[i].z += a1.z; }
		while( pos[i].y >= a2.y ){ pos[i].x -= a2.x; pos[i].y -= a2.y; pos[i].z -= a2.z; }
		while( pos[i].y < 0){ pos[i].x += a2.x; pos[i].y += a2.y; pos[i].z += a2.z; }
		while( pos[i].z >= a3.z ){ pos[i].x -= a3.x; pos[i].y -= a3.y; pos[i].z -= a3.z; }
		while( pos[i].z < 0){ pos[i].x += a3.x; pos[i].y += a3.y; pos[i].z += a3.z; }

		//wrong
		//while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a2.x; pos[i].z -= a3.x; }
		//while( pos[i].x < 0)     { pos[i].x += a1.x; pos[i].y += a2.x; pos[i].z += a3.x; }
		//while( pos[i].y >= a2.y ){ pos[i].x -= a1.y; pos[i].y -= a2.y; pos[i].z -= a3.y; }
		//while( pos[i].y < 0)     { pos[i].x += a1.y; pos[i].y += a2.y; pos[i].z += a3.y; }
		//while( pos[i].z >= a3.z ){ pos[i].x -= a1.z; pos[i].y -= a2.z; pos[i].z -= a3.z; }
		//while( pos[i].z < 0)     { pos[i].x += a1.z; pos[i].y += a2.z; pos[i].z += a3.z; }
	}
	return;
}

void Atoms::read_pos_5(ifstream &ifs, Vector3<double> &a1, Vector3<double> &a2, Vector3<double> &a3, double &lat_const)
{
	const double norm1 = INPUT.celldm1; 
	const double norm2 = INPUT.celldm2; 
	const double norm3 = INPUT.celldm3;
	assert(na>0);
	string tmp;
	pos = new Vector3<double>[na];
	posd = new Vector3<double>[na];
	allocate_pos = true;
	allocate_posd = true;
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

		// mohan added 2016-11-09
	//	while( pos[i].x >= norm1 ) pos[i].x -= norm1;	
	//	while( pos[i].x < 0 ) pos[i].x += norm1;
	//	while( pos[i].y >= norm2 ) pos[i].y -= norm2;	
	//	while( pos[i].y < 0 ) pos[i].y += norm2;
	//	while( pos[i].z >= norm3 ) pos[i].z -= norm3;	
	//	while( pos[i].z < 0 ) pos[i].z += norm3;
	
		// mohan updated 2018-04-26
		// correct
		while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a1.y; pos[i].z -= a1.z; }
		while( pos[i].x < 0){ pos[i].x += a1.x; pos[i].y += a1.y; pos[i].z += a1.z; }
		while( pos[i].y >= a2.y ){ pos[i].x -= a2.x; pos[i].y -= a2.y; pos[i].z -= a2.z; }
		while( pos[i].y < 0){ pos[i].x += a2.x; pos[i].y += a2.y; pos[i].z += a2.z; }
		while( pos[i].z >= a3.z ){ pos[i].x -= a3.x; pos[i].y -= a3.y; pos[i].z -= a3.z; }
		while( pos[i].z < 0){ pos[i].x += a3.x; pos[i].y += a3.y; pos[i].z += a3.z; }

		//wrong
		//while( pos[i].x >= a1.x ){ pos[i].x -= a1.x; pos[i].y -= a2.x; pos[i].z -= a3.x; }
		//while( pos[i].x < 0)     { pos[i].x += a1.x; pos[i].y += a2.x; pos[i].z += a3.x; }
		//while( pos[i].y >= a2.y ){ pos[i].x -= a1.y; pos[i].y -= a2.y; pos[i].z -= a3.y; }
		//while( pos[i].y < 0)     { pos[i].x += a1.y; pos[i].y += a2.y; pos[i].z += a3.y; }
		//while( pos[i].z >= a3.z ){ pos[i].x -= a1.z; pos[i].y -= a2.z; pos[i].z -= a3.z; }
		//while( pos[i].z < 0)     { pos[i].x += a1.z; pos[i].y += a2.z; pos[i].z += a3.z; }
	}
	return;
}

void Atoms::read_vel(ifstream &ifs)
{
	assert(na>0);
	this->vel = new Vector3<double>[na];
	allocate_vel = true;
//	cout << " atom number is " << na << endl;
	for(int i=0; i<na; ++i)
	{
		ifs >> vel[i].x >> vel[i].y >> vel[i].z;
		//cout << " vel=" << vel[i].x << " " << vel[i].y << " " << vel[i].z << endl;
	}
	return;
}


