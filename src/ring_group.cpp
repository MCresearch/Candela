#include "input.h"
#include "ring_group.h"

int RingGroup::nr_max = 50000;

RingGroup::RingGroup(){
	nr=0;
	na_of_ring=0;
}

RingGroup::~RingGroup()
{
}


bool RingGroup::push(short int *index, Atoms* atom)
{
	// compare to each number in the ring group
	for(int ir=0; ir<nr; ++ir) // number of rings stored in the RingGroup
	{
		if( IndexA_equals_RingB(index, r[ir]) ) 
		{
			return false;
		}
	}
	
	assert(nr<nr_max);
	assert(na_of_ring>0);

	this->r[nr].index_of_atoms = new int[this->na_of_ring];
	this->r[nr].na_ring = this->na_of_ring;
	for(int ia=0; ia<na_of_ring; ++ia)
	{
		r[nr].index_of_atoms[ia]=index[ia];
	}

	// print out
/*
	ofs_running << "Ring" << na_of_ring << endl;
	ofs_running << "1" << endl;
	ofs_running << INPUT.celldm1 * 0.529177 << " 0 0" << endl;
	ofs_running << "0 " << INPUT.celldm2 * 0.529177 << " 0" << endl;
	ofs_running << "0 0 " << INPUT.celldm3 * 0.529177 << endl;
	ofs_running << "O" << endl;
	ofs_running << na_of_ring << endl;
	ofs_running << "Cartesian" << endl;
	for(int i=0; i<na_of_ring; ++i)
	{
	double x=atom[index[i]].cx;
	double y=atom[index[i]].cy;
	double z=atom[index[i]].cz;
	while(x > INPUT.celldm1) x-=INPUT.celldm1;
	while(x < 0) x+=INPUT.celldm1;
	while(y > INPUT.celldm2) y-=INPUT.celldm2;
	while(y < 0) y+=INPUT.celldm2;
	while(z > INPUT.celldm3) z-=INPUT.celldm3;
	while(z < 0) z+=INPUT.celldm3;

	ofs_running  
	<< " " << x*0.529177 
	<< " " << y*0.529177 
	<< " " << z*0.529177 << endl;
	}
	// end print out
*/

	++nr;
	return true;
}

bool RingGroup::IndexA_equals_RingB(short int *indexA, const Ring &B)
{
	for(int ia=0; ia < na_of_ring; ++ia)
	{
		bool find=false;
		for(int ib=0; ib < na_of_ring; ++ib)
		{
			if( indexA[ia] == B.index_of_atoms[ib] )
			{
				find=true;
				break;
			}
		}
	
		if(!find)
		{
			return false;
		}
	}

	return true;
}

bool RingGroup::newbond(short int *indexA, const RingGroup *GB1, const RingGroup *GB2, const int &ngroup)
{
	assert(na_of_ring>0);

	if(INPUT.ring_definition==1) // Lixin's definition
	{
		bool* bond_exists = new bool[this->na_of_ring];
		for(int ia=0; ia<this->na_of_ring; ++ia)
		{
			bond_exists[ia] = false;
		}

		for(int ig=0; ig<ngroup; ++ig)
		{
			for(int ir=0; ir<GB1[ig].nr; ++ir)
			{
				IndexA_in_RingB(indexA, bond_exists, GB1[ig].r[ir]);
			}
		}

		for(int ia=0; ia<this->na_of_ring; ++ia)
		{
			if(bond_exists[ia] == false) 
			{
				delete[] bond_exists;
				return true;
			}
		}

		delete[] bond_exists;
		return false;
	}
	else if(INPUT.ring_definition==2) // Lingzhu's definition
	{
		for(int ig=0; ig<ngroup; ++ig)
		{	
			for(int ir=0; ir<GB2[ig].nr; ++ir)
			{
				if(IndexA_shortcutted_by_RingB(indexA, GB2[ig].r[ir]))
				{
					return false;
				}
			}
		}
		return true;
	}
}

bool RingGroup::IndexA_shortcutted_by_RingB(short int *indexA, const Ring &B)
{
	const int nover = count_overlap(indexA,B);
	if(na_of_ring==4 or na_of_ring==5)
	{
		// category 1
		if(B.na_ring==3 and nover==3 ) return true;
	}
	else if(na_of_ring==6 or na_of_ring==7)
	{
		// category 1
		if(B.na_ring==3 and nover==3 ) return true;
		if(B.na_ring==4 and nover==4 ) return true;
		// category 2
		if(B.na_ring==5 and nover==4 ) return true;
	}
	else if(na_of_ring==8 or na_of_ring==9)
	{
		// category 1
		if(B.na_ring==3 and nover==3 ) return true;
		if(B.na_ring==4 and nover==4 ) return true;
		if(B.na_ring==5 and nover==5 ) return true;
		// category 2
		if(B.na_ring==5 and nover==4 ) return true;
		if(B.na_ring==6 and nover==5 ) return true;
		// category 3
		if(B.na_ring==7 and nover==5 ) return true;
	}
	else if(na_of_ring==10 or na_of_ring==11)
	{
		// category 1
		if(B.na_ring==3 and nover==3 ) return true;
		if(B.na_ring==4 and nover==4 ) return true;
		if(B.na_ring==5 and nover==5 ) return true;
		if(B.na_ring==6 and nover==6 ) return true;
		// category 2
		if(B.na_ring==5 and nover==4 ) return true;
		if(B.na_ring==6 and nover==5 ) return true;
		if(B.na_ring==7 and nover==6 ) return true;
		// category 3
		if(B.na_ring==7 and nover==5 ) return true;
		if(B.na_ring==8 and nover==6 ) return true;
		// category 4
		if(B.na_ring==9 and nover==6 ) return true;
	}
	else if(na_of_ring==12 or na_of_ring==13)
	{
		// category 1
		if(B.na_ring==3 and nover==3 ) return true;
		if(B.na_ring==4 and nover==4 ) return true;
		if(B.na_ring==5 and nover==5 ) return true;
		if(B.na_ring==6 and nover==6 ) return true;
		if(B.na_ring==7 and nover==7 ) return true;
		// category 2
		if(B.na_ring==5 and nover==4 ) return true;
		if(B.na_ring==6 and nover==5 ) return true;
		if(B.na_ring==7 and nover==6 ) return true;
		if(B.na_ring==8 and nover==7 ) return true;
		// category 3
		if(B.na_ring==7 and nover==5 ) return true;
		if(B.na_ring==8 and nover==6 ) return true;
		if(B.na_ring==9 and nover==7 ) return true;
		// category 4
		if(B.na_ring==9 and nover==6 ) return true;
		if(B.na_ring==10 and nover==7 ) return true;
		// category 5
		if(B.na_ring==11 and nover==7 ) return true;
	}


	return false;
}

int RingGroup::count_overlap(short int *indexA, const Ring &B)
{
//	ofs_running << "new ring" << B.na_ring << " in ring" << na_of_ring << endl;
	int count=0;
	for(int ib=0; ib < B.na_ring; ++ib)
	{
		bool find=false;
		for(int ia=0; ia < this->na_of_ring; ++ia)
		{
			if(B.index_of_atoms[ib] == indexA[ia]) 
			{	
	//			ofs_running << "check " << ib << " " << ia << endl;
				find=true;
				break;
			}
		}
		if(find) ++count;
	}
	return count;
}

// 0519 new
void RingGroup::IndexA_in_RingB(short int *indexA, bool *bond_exists, const Ring &B)
{
	for(short int ia=0; ia < this->na_of_ring; ++ia)
	{
		bool find0=false;
		bool find1=false;
		bool find2=false;
		short int i1 = ia; 
		short int i2 = ia+1; if(i2>=na_of_ring) i2-=na_of_ring;
		short int i3 = ia+2; if(i3>=na_of_ring) i3-=na_of_ring;

		for(short int ib=0; ib < B.na_ring; ++ib)
		{
//-----------
//OPTION1 (not sure is correct)
//-----------
//			if( B.index_of_atoms[ib] == indexA[i1] ) find0=true;
//			else if( B.index_of_atoms[ib] == indexA[i2] ) find1=true;
//			else if( B.index_of_atoms[ib] == indexA[i3] ) find2=true;

//-----------
//OPTION2
//-----------
			short int j1 = ib; 
			short int j2 = ib+1; if(j2>=B.na_ring) j2-=B.na_ring;
			short int j3 = ib+2; if(j3>=B.na_ring) j3-=B.na_ring;

			if( B.index_of_atoms[j1] == indexA[i1] and
				B.index_of_atoms[j2] == indexA[i2] and
				B.index_of_atoms[j3] == indexA[i3])
			{
				bond_exists[ia]=true;
			}
			else if( B.index_of_atoms[j1] == indexA[i3] and
				B.index_of_atoms[j2] == indexA[i2] and
				B.index_of_atoms[j3] == indexA[i1])
			{
				bond_exists[ia]=true;
			}
// end of OPTION2

		}

		if(find0 == true and find1 == true and find2 == true) 
		{
			bond_exists[ia] = true;
		}
	}
	return;
}


// if ture IndexB is useless and should be get rid of. 
bool RingGroup::RingA_in_IndexB(const Ring &A, short int *indexB)
{

	for(int ia=0; ia < A.na_ring; ++ia)
	{
		bool find0=false;
		bool find1=false;
		bool find2=false;
		int i1 = ia; 
		int i2 = ia+1; if(i2>=A.na_ring) i2-=A.na_ring;
		int i3 = ia+2; if(i3>=A.na_ring) i3-=A.na_ring;

		for(int ib=0; ib < this->na_of_ring; ++ib)
		{
			if( A.index_of_atoms[i1] == indexB[ib] ) find0=true;
			if( A.index_of_atoms[i2] == indexB[ib] ) find1=true;
			if( A.index_of_atoms[i3] == indexB[ib] ) find2=true;
		}

		if(find0 == false and find1 == false and find2 == false) 
		{
			return false; // then we should count ringB
		}
	}

	return true;


/*
	for(int ia=0; ia < A.na_ring; ++ia)
	{
		bool find0=false;
		bool find1=false;
		bool find2=false;
		int i1 = ia; 
		int i2 = ia+1; if(i2>=A.na_ring) i2-=A.na_ring;
		int i3 = ia+2; if(i3>=A.na_ring) i3-=A.na_ring;

		for(int ib=0; ib < this->na_of_ring; ++ib)
		{
			if( A.index_of_atoms[i1] == indexB[ib] ) find0=true;
			if( A.index_of_atoms[i2] == indexB[ib] ) find1=true;
			if( A.index_of_atoms[i3] == indexB[ib] ) find2=true;
		}

		if(find0 == true and find1 == true and find2 == true) 
		{
			return true; // then we should not count ringB
		}
	}

	return false; 
*/
}

bool RingGroup::sumUp(Atoms* atom, short int *indexB)
{
	double sumx = 0.0;
	double dx;

	short int i1, i2;
	for(short int ia=0; ia < this->na_of_ring; ++ia)
	{
		i1=indexB[ia];
		if(ia!=na_of_ring-1) i2 = indexB[ia+1];
		else i2=indexB[0];
		dx = atom[i1].cx - atom[i2].cx;
		while( dx >  Atoms::celldm1_half ) dx -= Atoms::celldm1;
		while( dx < -Atoms::celldm1_half ) dx += Atoms::celldm1;
		sumx += dx;
	}
	
	if( abs(sumx) > 1.0e-6 ) return false;

	double sumy = 0.0;
	double dy;

	for(short int ia=0; ia < this->na_of_ring; ++ia)
	{
		i1=indexB[ia];
		if(ia!=na_of_ring-1) i2 = indexB[ia+1];
		else i2=indexB[0];
		dy = atom[i1].cy - atom[i2].cy;
		while( dy >  INPUT.celldm2/2.0 ) dy -= Atom::celldm2;
		while( dy < -INPUT.celldm2/2.0 ) dy += Atom::celldm2;
		sumy += dy;
	}

	if( abs(sumy) > 1.0e-6 ) return false;

	double sumz = 0.0;
	double dz;

	for(short int ia=0; ia < this->na_of_ring; ++ia)
	{
		i1=indexB[ia];
		if(ia!=na_of_ring-1) i2 = indexB[ia+1];
		else i2=indexB[0];
		dz = atom[i1].cz - atom[i2].cz;
		while( dz >  Atom::celldm3_half ) dz -= Atom::celldm3;
		while( dz < -Atom::celldm3_half ) dz += Atom::celldm3;
		sumz += dz;
	}

	if( abs(sumz) > 1.0e-6 ) return false;

//	cout << "sum " << sumx << " " << sumy << " " << sumz << endl;
	return true;
}
