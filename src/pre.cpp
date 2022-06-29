#include "cellFile.h"
#include "input.h"
#include "pre.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"

PRE::PRE()
{
}

PRE::~PRE(){}

void PRE::Routine()
{
	TITLE("PRE","Routine");
	
	cout << "Compute the Presolvation Structure of Water" << endl;

	//************
	// OPEN FILES 
	//************
	PT.setup_PT();

	// "nion": find out how many ions are involved in becoming an ion during the trajectory
	// the indexes are saved in "ion_pt"
	bool* already_found = new bool[PT.npt];
	for(int i=0; i<PT.npt; ++i) already_found[i]=false;
	this->nion=0;
	for(int i=0; i<PT.npt; ++i)
	{
		if(already_found[i]==true) continue;
		//cout << PT.ion_index[i] << endl;
		++nion;
		already_found[i]=true;
		for(int j=i+1; j<PT.npt; ++j)
		{
			if(PT.ion_index[j]==PT.ion_index[i])
			{
				already_found[j]=true;
			}
		}
	}
	cout << "number of different ions is " << nion << endl;
	delete[] already_found;

	// again, fill ion_pt
	ion_pt = new int[nion]();
	for(int ia=0; ia<nion; ++ia)
	{
		ion_pt[ia]=-1;
	}

	//----------------------------------------------------
	// ion_pt: index for all atoms that can be ions 
	// PT.npt: number of proton transfer events
	//----------------------------------------------------
	int k=0;
	for(int ip=0; ip<PT.npt; ++ip)
	{
		bool sameion = false;
		for(int ia=0; ia<k; ++ia)
		{
			if(ion_pt[ia] == PT.ion_index[ip])
			{
				sameion=true;
			}
		}
		if(!sameion)
		{
			ion_pt[k]=PT.ion_index[ip];
			cout << k << " " << ion_pt[k] << endl;
			++k;
		}
	}
	assert(k==nion);
	cout << "nion: " << nion << endl;

	// initialize
	int *ion_nacc = new int[nion]();
	int *ion_ndon = new int[nion]();

	//**************************
	// BEGIN CALCULATING DATA
	//**************************
	this->count_geometry_number=0;
	int ipt=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

//		if(cel.snapshot_time > 1.0) continue; //mohan test

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;
		search_presolvation(cel, igeo, ion_nacc, ion_ndon);
	}	

	//********************
	//	PRINT INFORMATION
	//********************
	delete[] ion_pt;
	delete[] ion_nacc;
	delete[] ion_ndon;

	return;
}


void PRE::search_presolvation(const Cell &cel, const int &igeo, int* ion_nacc, int* ion_ndon)
{
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


	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

	HBs::setup_water(cel, water);
	ofs_running << "snapshot " << cel.snapshot_index << " nions = " << Water::nions << endl;

	//cout << cel.snapshot_index << endl;

	//-------------------------------------------------
	// all pt events included, including rattling
	//-------------------------------------------------
	for(int ip=1; ip<PT.npt; ++ip)
	{
		if(cel.snapshot_index == PT.snapshot_index[ip])
		{
			if(PT.type_pt[ip-1]!="rattling")
			{
				ofs_running << "PT event " << setw(5) << ip+1 
				<< setw(5) << PT.ion_index[ip] 
				<< setw(12) << cel.snapshot_index 
				<< setw(12) << cel.snapshot_time;

				// find the solvation structure of the last ion (I0) 
				int find2=-1;
				for(int ia2=0; ia2<nion; ++ia2)
				{
					if(ion_pt[ia2]==PT.ion_index[ip-1]) 
					{
						find2=ia2;
					}
				}		
				assert(find2!=-1);
				// because ion_index starts from 1, so here -1 is needed.
				const int ia3 = PT.ion_index[ip-1]-1;


				// find the structures of ion (N1)
				int find=-1;
				for(int ia2=0; ia2<nion; ++ia2)
				{
					if(ion_pt[ia2]==PT.ion_index[ip]) 
					{
						find=ia2;
					}
				}		
				assert(find!=-1);

				// I0 is the structure index before this PT (ip-1), when the ion (currently) was still not an ion
				ofs_running << " I0 " << ion_nacc[find2] << ion_ndon[find2];
				ofs_running << " N0 "<< ion_nacc[find] << ion_ndon[find];
				ofs_running << " I1 " << water[ia3].naccept << water[ia3].ndonate; 
				// N1 is the structure index after this PT (ip), now it is the ion
				ofs_running << " N1 " << PT.nacc[ip]  << PT.ndon[ip];




				if(ion_nacc[find2]==PT.nacc[ip] and ion_ndon[find2]==PT.ndon[ip])
				{
					ofs_running << " PRE_YES";
				}
				else
				{
					ofs_running << " PRE_NO";
				}

				ofs_running << " " << PT.type_pt[ip-1] << " ";

//				ofs_running << " " << cel.atom[ito].pos[PT.ion_index[ip]-1].z;
				ofs_running << endl;
			}
		}
	}


	// only if there is only one ion in the water
	if(Water::nions==1)
	{
		// save solvation structures (accepted number of HBs and donating number 
		// of HBs) for all ions that appear in the whole trajectory
		for(int ia=0; ia<nion; ++ia)
		{
			int ia2 = ion_pt[ia]-1;
			assert(ia2>=0);
//			cout << ia2 << " " << water[ia2].naccept << endl;
			ion_nacc[ia] = water[ia2].naccept;
			ion_ndon[ia] = water[ia2].ndonate;
		}	
	}

	delete[] water;

}
