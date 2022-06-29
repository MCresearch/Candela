#include "eig.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Eig::Eig() 
{
}

Eig::~Eig() 
{
}

void Eig::Routine()
{
	cout << "===== Compute the distribution of Eigenvalues =====" << endl;	

	// plot out the final distance between H and wannier centers
	this->dr=INPUT.dr;
	assert(dr<0);
	double rcut=INPUT.rcut;
	assert(rcut<0);
	this->nr=abs(rcut/dr)+10; //+10 for safety
	assert(nr>0);
	this->dis_eig = new double[nr]();



	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		cout << "igeo=" << igeo << endl;

		// cel : input geometry file
		CellFile cel;

		//ofs_running << "geometry " << igeo%INPUT.geo_interval << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;


		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;

		read_eig(cel, dis_eig);
	}	

	// print out the distribution of eigenvalues
	ofstream ofs("dis_eig.dat");

	double sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_eig[ir]*dr;
	}
	//cout << "sum=" << sum << endl;

	if(sum!=0)
	{
		for(int ir=0; ir<nr; ++ir)
		{
			dis_eig[ir]/=sum;
		}
	}
	else
	{
		cout << " Error! sum=0.0" << endl;
		exit(0);
	}

	sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_eig[ir]*dr;
	}
	//cout << "sum is " << sum << endl;

	for(int ir=0; ir<nr; ++ir)
	{
		ofs << dr*(ir+1) << " " << -dis_eig[ir] << endl;
	}
	ofs.close();

	delete[] dis_eig;

	return;
}





void Eig::read_eig(Cell &cel, double* dist_eig)
{
//	cout << "reading Wannier centers" << endl;

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

	// after the shift, only those eigenvalues < 0 are recorded 
	for(int ib=0; ib<INPUT.nbands; ++ib)
	{
		cel.eig[ib] -= INPUT.shift;	
	}

	for(int ib=0; ib<INPUT.nbands; ++ib)
	{
		int index=cel.eig[ib]/this->dr;
		if(index<nr and index>=0) 
		{
			dis_eig[index]++;
		}
	}

	return;
}
