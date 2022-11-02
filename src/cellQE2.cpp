#include "cellFile.h"
#include "input.h"
#include "gfun.h"
#include "const.h"
//qianrui
void lim2cel(double &x);
bool CellFile::ReadGeometry_QE2( Cell &cel, ifstream &ifs )
{
	TITLE("CellFile","ReadGeometry_QE2");
	//celldm
	if(CellFile::first_read)
	{
		string useless;
		while(ifs>>useless)
		{
			if(useless=="celldm(1)=")
			{
				ifs>>CellFile::celldm;
				CellFile::celldm *= P_BOHR;
				if(RANK == 0) cout<<"celldm1: "<<CellFile::celldm
				                  <<"; celldm2: "<<CellFile::celldm
								  <<"; celldm3: "<<CellFile::celldm<<endl;
				break;
			}
		}
		if(ifs.eof())
		{
			cout<<"Geo isn't enough!"<<endl;
			exit(0);
		}
		CellFile::first_read = false;
	}
	const int ntype = INPUT.ntype;
	double celldm1, celldm2, celldm3;
	if(INPUT.celldm1 * INPUT.celldm2 * INPUT.celldm3 > 1e-4)
	{
		celldm1 =INPUT.celldm1;
		celldm2 =INPUT.celldm2;
		celldm3 =INPUT.celldm3;
	}
	else
	{
		celldm1 = celldm2 = celldm3 = CellFile::celldm;
	}
	
	cel.nat=INPUT.natom;
	cel.a1.x=celldm1;	cel.a1.y=0;		cel.a1.z=0;
	cel.a2.x=0;		cel.a2.y=celldm2;	cel.a2.z=0;
	cel.a3.x=0;		cel.a3.y=0;		cel.a3.z=celldm3;
	cel.volume=celldm1*celldm2*celldm3;
	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	cout << "count geometry " << count_geometry << endl;
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt; 
	delete[] cel.atom;
	cel.atom = new Atoms[ntype];

	string useless;
	string txt;
	// locate band info
	if(INPUT.nbands>0)
	{
		while(getline(ifs,txt))
		{
			if(txt.find("bands (ev):")!=string::npos)
			{
			//	cout << txt << endl;
				getline(ifs,txt);
			//	cout << txt << endl;
				break;
			}
		}
		cel.nbands = INPUT.nbands;
		delete[] cel.eig;
		cel.eig = new double[INPUT.nbands]();
		//ifs >> useless;
		for(int iband=0; iband<INPUT.nbands; iband++)
		{
			ifs >> cel.eig[iband];
		}	
	}
	
	//locate atomic positions
	while(ifs>>useless)
	{
		getline(ifs,txt);
		if(useless=="ATOMIC_POSITIONS")
		{
			if(txt == " (crystal)") INPUT.cartesian = false;
			else					INPUT.cartesian = true;
			break;
		}
	}
	if(ifs.eof())
	{
		cout<<"Geo isn't enough!"<<endl;
		return false;
	}
	if(ntype==1)
	{
		cel.atom[0].na = cel.nat;
		cel.atom[0].id = INPUT.id1;
	}
	else if(ntype==2)
	{
		cel.atom[0].na = INPUT.natom1;
		cel.atom[1].na = INPUT.natom2;
		cel.atom[0].id = INPUT.id1;
		cel.atom[1].id = INPUT.id2;
	}
	else if(ntype==3)
	{
		cel.atom[0].na = INPUT.natom1;
		cel.atom[1].na = INPUT.natom2;
		cel.atom[2].na = INPUT.natom3;
		cel.atom[0].id = INPUT.id1;
		cel.atom[1].id = INPUT.id2;
		cel.atom[2].id = INPUT.id3;
	}


	// mohan updated 2018-12-14
	for(int it=0; it<ntype; ++it)
	{
		if(cel.atom[it].id=="O")
		{
			cel.atom[it].mass=15.9994;
		}
		else if(cel.atom[it].id=="H")
		{
			cel.atom[it].mass=1.00794;
		}
		else if(cel.atom[it].id=="D")
		{
			cel.atom[it].mass=2.014;
		}
		else if(cel.atom[it].id=="Al")
		{
			cel.atom[it].mass=26.981539;
		}
		else if(cel.atom[it].id=="Mg")
		{
			cel.atom[it].mass=24.305;
		}
		else if(cel.atom[it].id=="Be")
		{
			cel.atom[it].mass=9.012182;
		}
	}



	double tmpx,tmpy,tmpz;
	for(int it=0; it<ntype; ++it)
        {
                delete[] cel.atom[it].pos;
        }// renxi added 20200902
	int* record_na = new int[ntype];
	for(int it=0; it<ntype; ++it)
	{
		record_na[it] = 0;
		cel.atom[it].pos = new Vector3<double>[cel.atom[it].na];
	}
	for(int ia=0; ia<INPUT.natom; ++ia)
	{
		ifs>>useless>>tmpx>>tmpy>>tmpz;
		for(int it=0; it<ntype; ++it) // renxi fixed 20200902
		{
			if (useless == cel.atom[it].id || ntype == 1)
			{	
				int ia1 = record_na[it];
  				lim2cel(tmpx);
				lim2cel(tmpy);
				lim2cel(tmpz);
				if(INPUT.cartesian)
                {
					//cout << tmpx << " " << tmpy << " " << tmpz << endl;
					cel.atom[it].pos[ia1].x=tmpx;
					cel.atom[it].pos[ia1].y=tmpy;
					cel.atom[it].pos[ia1].z=tmpz;
                }
				else
				{
					cel.atom[it].pos[ia1].x=tmpx*celldm1;
					cel.atom[it].pos[ia1].y=tmpy*celldm2;
					cel.atom[it].pos[ia1].z=tmpz*celldm3;
				}
				//cout << cel.atom[it].pos[ia1].x << " " << cel.atom[it].pos[ia1].y << " " << cel.atom[it].pos[ia1].z << endl;
				record_na[it]++;
				//cout<<cel.atom[it].pos[ia].x<<' '<<cel.atom[it].pos[ia].y<<' '<<cel.atom[it].pos[ia].z<<endl;
			}
		}
	}
	delete[] record_na;
	

//	cout << 1 << endl;
	return true;
}

void lim2cel(double &x)
{
	double celldm1 = INPUT.celldm1;
	if ((INPUT.celldm1 != INPUT.celldm2 or INPUT.celldm1 != INPUT.celldm3) and INPUT.length_unit == "angstrom")
	{
		cout << "The function lim2cel in cellQE2.cpp requires all celldm to be equal." << endl;
		exit(0);
	}
	if(INPUT.cartesian)
	{
		while(x>=celldm1) x-=celldm1;
		while(x<0) x+=celldm1;
	}
	else
	{
		while(x>=1) x-=1;
		while(x<0) x+=1;
	}
}


