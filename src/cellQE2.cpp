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
	cel.cal_volume();
	static int count_geometry=0;
	cel.snapshot_index = count_geometry;
	++count_geometry;
	cout << "count geometry " << count_geometry << endl;
	cel.snapshot_time = cel.snapshot_index * INPUT.msd_dt;

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


	//read atom positions
	double tmpx,tmpy,tmpz;
	int* record_na = new int[ntype];
	ZEROS(record_na, ntype);
	for(int ia=0; ia<INPUT.natom; ++ia)
	{
		ifs>>useless>>tmpx>>tmpy>>tmpz;

		//get which it
		int it = -1;
		if(ntype == 1)
		{
			it = 0;
			cel.atom[0].id = useless;
		}
		else
		{
			for(int it0=0; it0<ntype; ++it0)
			{
				if (useless == cel.atom[it0].id)
				{
					it = it0;
				}
			}
		}
		if(it < 0)
		{
			cout<<"We do not find atom id "<<useless<<endl;
			exit(0);
		}

		//get positions
		int ia1 = record_na[it];
  		lim2cel(tmpx);
		lim2cel(tmpy);
		lim2cel(tmpz);
		if(INPUT.cartesian)
        {
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
		record_na[it]++;
	}
	for(int it=0; it<ntype; ++it)
	{
		assert(record_na[it] == cel.atom[it].na);
	}
	delete[] record_na;
	cel.atom_mass();

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


