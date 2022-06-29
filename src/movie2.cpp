#include "movie2.h"
#include "input.h"
#include "HBs.h"

void movie2::Routine()
{	
	ofs_running << "Produce a movie for MD trajectory" << endl;

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;
	ifstream ifs("Oindex.txt");
	// only for double jump, the first index is the possible accepted water's O.
	// The second index is hydronium, the third and fourth are destinations respectively for double jumps.
	int* Ochain = new int[INPUT.natom1]();

	for(int ia=0; ia<INPUT.u1; ia++)
	{
		ifs >> Ochain[ia];
		Ochain[ia]--;
		cout << ia << " " << Ochain[ia] << endl;
	}
	ifs.close();
	
	// output data
	ofstream ofs(INPUT.geo_out.c_str());
	bool init = false;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//ofs_running << "geometry " << igeo%INPUT.geo_interval << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
       		 {
           		 cel.read_and_used=false;
       		 }
		else cel.read_and_used=true;

		cout << igeo << " " << cel.read_and_used << endl;


		stringstream ss; 
		ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		ofs_running << "igeo=" << igeo << endl;
		cout << "igeo=" << igeo << endl;
		static int* Hindex = new int[INPUT.natom2]();
		static int nO;
		static int nH;
		snapshot(ofs,cel,igeo, Ochain, Hindex, init, nO, nH);
	}

	ofs.close();


	return;
}

void movie2::snapshot(ofstream &ofs,Cell &cel,int &igeo, int* &Ochain, int* &Hindex, bool &init, int &nO, int &nH)
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

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];
	HBs::setup_water(cel, water);

	
	if(not init)
	{
		nO = INPUT.u1;
		bool already_in = false;
	/*	for(int idon=0; idon<water[Ochain[1]].ndonate; idon++)
		{
			int index = water[Ochain[1]].donateO[idon];
			for(int ia=0; ia<nO; ia++)
			{
				if(Ochain[ia] == index)
				{
					already_in = true;
					break;
				}
			}
			if(already_in)
			{
				already_in = false;
				continue;
			}
			Ochain[nO] = index;
			cout << nO << " " << index << endl;
			nO++;
		}
*/
		for(int isa = 0; isa < INPUT.u1; isa++)
		{
			already_in = false;
			for(int idon=0; idon<water[Ochain[isa]].ndonate; idon++)
			{
				int index = water[Ochain[isa]].donateO[idon];
				for(int ia=0; ia<nO; ia++)
				{
					if(Ochain[ia] == index)
					{
						already_in = true;
						break;
					}
				}
				if(already_in)
				{
					already_in = false;
					continue;
				}
				Ochain[nO] = index;
				cout << nO << " " << index << endl;
				nO++;
			}
			already_in = false;
			for(int idon=0; idon<water[Ochain[isa]].naccept; idon++)
			{
				int index = water[Ochain[isa]].acceptO[idon];
				for(int ia=0; ia<nO; ia++)
				{
					if(Ochain[ia] == index)
					{
						already_in = true;
						break;
					}
				}
				if(already_in)
				{
					already_in = false;
					continue;
				}
				Ochain[nO] = index;
				cout << nO << " " << index << endl;
				nO++;
			}
		}
		ofstream ofs_totalOlist("totalOindex.txt");
		for(int ia=0; ia<nO; ia++)
		{
			cout << ia << " " << Ochain[ia]	<< endl;
			ofs_totalOlist << ia << " " << Ochain[ia] << endl;
		}
		ofs_totalOlist.close();
		int iH = 0;
		already_in = false;

		for(int ia=0; ia<nO; ia++)
		{
			int index = Ochain[ia];
			cout << water[index].nH << endl;
			for(int iah = 0; iah<water[index].nH; iah++)
			{
				int newH = water[index].indexH[iah];
				cout << index << " " << newH << endl;
				for(int iih=0; iih<iH; iih++)// check whether newH is already documented
				{
					if(Hindex[iih] == newH)
					{
						already_in = true;
						break;
					}
				}
				if(already_in)
				{
					already_in = false;
					continue;
				}
				Hindex[iH] = newH;
				iH++;
			}
		}
		assert(iH == nO*2+1);
		nH = iH;
		init = true;
	}// end init
	for(int ia=0; ia<nH; ia++)
	{
		cout << Hindex[ia] << " ";
	}
	cout << endl;
	ofs << nO+nH << endl;
	ofs << cel.snapshot_index << " " << cel.snapshot_time << endl;
	for(int ia=0; ia<nO; ia++)
	{
		ofs << "O " << cel.atom[ito].pos[Ochain[ia]].x << " " << cel.atom[ito].pos[Ochain[ia]].y << " " << cel.atom[ito].pos[Ochain[ia]].z << endl;
	}
	for(int ia=0; ia<nH; ia++)
	{
		ofs << "H " << cel.atom[ith].pos[Hindex[ia]].x << " " << cel.atom[ith].pos[Hindex[ia]].y << " " << cel.atom[ith].pos[Hindex[ia]].z << endl;
	}
		
	return;
}

