#ifndef RING_GROUP_H
#define RING_GROUP_H

#include "gfun.h"
#include "atoms.h"

// define a ring here
class Ring
{
	public:
    
	Ring(){};
    ~Ring(){};
	
	int na_ring;
	int* index_of_atoms;


};

// define a group of rings here
class RingGroup
{
	public:
    
	RingGroup();
    ~RingGroup();

	void Init();
	Ring* r;
	int nr; // number of rings
	int na_of_ring;
	static int nr_max; 

	bool push(short int* index, Atoms* atom);

	bool newbond(short int *indexA, const RingGroup *GB1, const RingGroup *GB2, const int &ngroup);
	void IndexA_in_RingB(short int *indexA, bool *bond_exists, const Ring &B);

	bool sumUp(Atoms* atom, short int *indexB);


	private:

	bool IndexA_equals_RingB(short int *indexA, const Ring &B);
	bool RingA_in_IndexB(const Ring &A, short int *indexB);
	bool IndexA_shortcutted_by_RingB(short int *indexA, const Ring &B);

	// method by Lingzhu
	int count_overlap(short int *indexA, const Ring &B);

};

#endif //RINGGROUP
