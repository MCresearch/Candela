#include "input.h"
#include "water.h"
#include "math.h"

int Water::nions=0;

Water::Water()
{
	indexO=-1;
	indexH = new int[INPUT.nHB_max];
	disH = new double[INPUT.nHB_max];
	acceptO = new int[INPUT.nHB_max];
	acceptH = new int[INPUT.nHB_max];
	donateO = new int[INPUT.nHB_max];
	donateH = new int[INPUT.nHB_max];
	accept_angle = new double[INPUT.nHB_max];
	accept_disO = new double[INPUT.nHB_max];
	donate_angle = new double[INPUT.nHB_max];
	donate_disO = new double[INPUT.nHB_max];
	for(int i=0; i<INPUT.nHB_max; ++i)
	{
		indexH[i]=-1;	
		disH[i]=-1.0;
		acceptO[i]=-1;
		acceptH[i]=-1;
		donateO[i]=-1;
		donateH[i]=-1;
		accept_angle[i]=-1;
		donate_angle[i]=-1;
		accept_disO[i]=-1;
		donate_disO[i]=-1;
	}
	nH = 0;
	nO = 0;
	nO_H = 0;
	naccept = 0;
	ndonate = 0;
	for(int i=0; i<3; ++i)
	{
		dipole[i]=0.0;
	}
	dipole_sum = 0.0;
}

Water::~Water()
{

	delete[] indexH;
	delete[] disH;
	delete[] acceptO;
	delete[] acceptH;
	delete[] donateO;
	delete[] donateH;
	delete[] accept_angle;
	delete[] donate_angle;
	delete[] accept_disO;
	delete[] donate_disO;

}
