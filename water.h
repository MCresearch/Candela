#ifndef WATER_H
#define WATER_H


class Water
{
	public:

	Water();
	~Water();

	static int nions; // number of ions, added 2016-11-28
	int indexO; // starting from 0
	int nH; // number of adjacent H atoms 
	int nO; // number of adjacent O atoms
	int nO_H; // number of adjacent atoms for H
	int naccept;
	int ndonate; 

	int* indexH; // index of H atoms, starting from 0
	double* disH; // in Angstroms

	int* acceptO; // indexes of accepted O
	int* acceptH; // indexes of accepted H
	double* accept_angle; // angles H-O-O 
	double* accept_disO;

	int* donateO; // indexes of donating O
	int* donateH; // indexes of donating O
	double* donate_angle; // angles H-O-O
	double* donate_disO; 

	double bond_angle; // bond angle of this water

	double dipole[3]; // dipole moment of this water, mohan added 2017-02-28
	double dipole_sum; // dipole moment of this water, mohan added 2017-02-28

	bool used;
	


	private:

};

#endif
