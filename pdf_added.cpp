#include "cellFile.h"
#include "input.h"
#include "pdf_added.h"
#include "math.h"
#include "HBs.h"
#include "Honeycutt.h"

bool PDF_ADDED::atom_in_ion2(const Cell &cel, const int &it1, const int &ia, const int &it2, const int &ia2)
{
	if(cel.atom[it1].id=="H" or cel.atom[it1].id=="D")
	{
		double ho_dis = distance(cel.atom[it1].pos[ia],cel.atom[it2].pos[ia2],
			cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
		if(ho_dis<INPUT.rcut_oh)
		{
			return false; // shoult not use this oxygen (it2, ia2)
		}
		else
		{
			return true;
		}
	}

	return true;
}


bool PDF_ADDED::atom_in_ion(const Cell &cel, const int &it_in, const int &ia, int &io_of_ion)
{
	int ito=-1;
	int ith=-1;
	int itc=-1;
	int itcl=-1;
	int itna=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
		else if(cel.atom[it].id=="Na") itna=it;
	}


	// to decide whether O belongs to hydronium(H3O+) or hydroxide(OH-)
	if(cel.atom[it_in].id=="O")
	{
		// this is the index of ion
		io_of_ion = ia;
		
		int nH=0;
		for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
		{
			double ho_dis = distance(cel.atom[ito].pos[ia],cel.atom[ith].pos[ia2],
					cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
			if(ho_dis<INPUT.rcut_oh)
			{
				++nH;
			}
		}

		if(INPUT.system=="hydroxide") 
		{
			if(nH==1) return true;
			else return false;
		}
		else if(INPUT.system=="hydronium")
		{
			if(nH==3) return true;
			else return false;
		}
	}
	else if(cel.atom[it_in].id=="H" or cel.atom[it_in].id=="D")
	{
		int belong_to = -1;
		for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
		{
			double ho_dis = distance(cel.atom[ith].pos[ia],cel.atom[ito].pos[ia2],
					cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
			if(ho_dis<INPUT.rcut_oh)
			{
				belong_to=ia2;
				break;
			}
		}


		if(belong_to==-1) return false;
		else io_of_ion = belong_to;
		
		if(belong_to>=0)
		{
			int nH=0;
			for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
			{
				double ho_dis = distance(cel.atom[ito].pos[belong_to],cel.atom[ith].pos[ia2],
						cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
				if(ho_dis<INPUT.rcut_oh)
				{
					++nH;
				}
			}

			if(INPUT.system=="hydroxide") 
			{
				if(nH==1) return true;
				else return false;
			}
			else if(INPUT.system=="hydronium")
			{
				if(nH==3) return true;
				else return false;
			}

		}

	}

}


// compute 'delta', which is the proton transfer coordinates
double PDF_ADDED::compute_delta(const Cell &cel, const Water* water, const int &it_in, const int &ia)
{
	// elements
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

	// the return value
	double delta=10000.00;

	if(cel.atom[it_in].id=="H" or cel.atom[it_in].id=="D")
	{
		int belong_to = -1;
		for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
		{
			double ho_dis = distance(cel.atom[ith].pos[ia],cel.atom[ito].pos[ia2],
						cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
			if(ho_dis<INPUT.rcut_oh)
			{
				belong_to=ia2;
				break;
			}
		}

		if(belong_to==-1) return false;
		
		if(belong_to>=0)
		{
			// adjacent oxygen atoms of hydronium or hydroxide
			for(int io=0; io<water[belong_to].naccept; ++io)
			{
				const int oindex = water[belong_to].acceptO[io];
				// attached hydrogen atoms of adjacent water molecules
				for(int ih=0; ih<water[oindex].nH; ++ih)
				{
					const int hindex = water[oindex].indexH[ih];	
					double hox_dis = distance(cel.atom[ith].pos[hindex],cel.atom[ito].pos[belong_to],
									cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
					double ho_dis = distance(cel.atom[ith].pos[hindex],cel.atom[ito].pos[oindex],
									cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
					double tmp = abs(hox_dis-ho_dis);
					if(tmp<delta)
					{
						delta = tmp;
					}
				}
			}
		}
	}
	else if(cel.atom[it_in].id=="O")
	{
		int belong_to=ia;
		// index of adjacent oxygen atoms of hydronium or hydroxide
		for(int io=0; io<water[belong_to].naccept; ++io)
		{
			// index of adjacent atoms of water ions
			const int oindex = water[belong_to].acceptO[io];
			// indexed of hydrogen atoms from adjacent water molecules
			for(int ih=0; ih<water[oindex].nH; ++ih)
			{
				const int hindex = water[oindex].indexH[ih];	
				// distance of H from adjacent O to ion
				double hox_dis = distance(cel.atom[ith].pos[hindex],cel.atom[ito].pos[belong_to],
						cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
				// distance of H from its host O
				double ho_dis = distance(cel.atom[ith].pos[hindex],cel.atom[ito].pos[oindex],
						cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
				// if tmp < 0, the H is more close to ion
				// if tmp > 0, the H is more close to its host O
				double tmp = abs(hox_dis-ho_dis);
				if(tmp<delta)
				{
					delta = tmp;
				}
			}
		}
	}


	return delta;
}



// return value: true means skip 
bool PDF_ADDED::water_ions(const Cell &cel, const Water* water, const int &io_of_ion, const int &ia2, bool &should_count)
{
	if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		if(INPUT.func==1)
		{
			// criterion 1
			// to exclude O* from H*O
			// but for O*O, this should not be performed
			// if(!PDF_ADDED::atom_in_ion2(cel,it,ia,it2,ia2)) continue;
			// criterion 2
			// to separate HBs (3,4,or more) from the whole g(r)
			if(INPUT.nacc>=0)
			{
				if(INPUT.nacc < 10)
				{
					if(water[io_of_ion].naccept!=INPUT.nacc) return 1;
				}
				else if(INPUT.nacc >= 10) // for example, 40 means >=4
				{
					if(water[io_of_ion].naccept < INPUT.nacc/10) return 1;
				}
			}
		}
		else if(INPUT.func==2) // mohan add 2016-12-01
		{
			if(INPUT.nacc>=0)
			{
				if(INPUT.nacc < 10)
				{
					if(water[io_of_ion].naccept!=INPUT.nacc) return 1;
				}
				else if(INPUT.nacc >= 10) // for example, 40 means >=4
				{
					if(water[io_of_ion].naccept < INPUT.nacc/10) return 1;
				}
			}
			bool one_of_acceptedO=false;
			for(int ia3=0; ia3<water[io_of_ion].naccept; ++ia3)
			{
				if(ia2==water[io_of_ion].acceptO[ia3])
				{
					one_of_acceptedO=true;
					break;
				}
			}
			if(one_of_acceptedO==false) return 1;
		}
		// mohan added 2017-02-19
		else if(INPUT.func==3) // only those O*(O_Donate_to_O*) recorded
		{
			bool find_donate=false;
			for(int idd=0; idd<water[io_of_ion].ndonate; ++idd)
			{
				if(water[io_of_ion].donateO[idd]==ia2)
				{
					find_donate=true;
				}
			}
			if(find_donate) return 1;
		}
		should_count=true;
	}// end system hydronium and hydroxide

	return 0;
}
