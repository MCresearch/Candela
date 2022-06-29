#include "cellFile.h"
#include "input.h"
#include "mdp.h"
#include "math.h"
#include "water.h"
#include "HBs.h"

MDP::MDP()
{
	ave_donate = 0.0;
	ave_accept = 0.0;
	ave_count = 0;
}

MDP::~MDP(){}

void MDP::Routine()
{
	TITLE("MDP","Routine");
	
	ofs_running << "Compute Mean Density Profile based on Instantaneous Liquid Interface (ILI)" << endl;
	ofs_running << "Recommend paper: J. Phys. Chem. B 2010, 114, 1954-1958." << endl;

	assert(INPUT.dr>0.0);
	assert(INPUT.rcut>0.0);

	this->dr = INPUT.dr;
	this->rcut = INPUT.rcut;
	this->nmesh = int(rcut / dr) +  1;

	// record the information for ion
	this->record_ion = new double*[7];
	for(int i=0; i<7; ++i) record_ion[i] = new double[nmesh]();

	// record the information for ion
	this->record_water = new double*[10];
	for(int i=0; i<10; ++i) record_water[i] = new double[nmesh]();


	// allocate for joint_cd (joint conditional distribution)
	assert(INPUT.u1>0);
	assert(INPUT.u2>0);
	this->joint_cd = new double*[INPUT.u1];
	for(int i=0; i<INPUT.u1; ++i)
	{
		joint_cd[i] = new double[INPUT.u2]();
	}
	assert(INPUT.u1>0);
	assert(INPUT.u2>0);
	this->du1 = 2.0/INPUT.u1;
	this->du2 = 2.0/INPUT.u2;
	cout << "du1=" << du1 << endl;
	cout << "du2=" << du2 << endl;
	

	// setup interface for each snapshot and its
	// associated gradient
	this->interface = new double*[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->interface[ix] = new double[INPUT.ny]();
	}
	this->gradient = new double**[INPUT.nx];
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		this->gradient[ix] = new double*[INPUT.ny];
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			this->gradient[ix][iy] = new double[3]();
		}
	}

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	// input ili file
	ifstream ifs(INPUT.ili_file.c_str());
	if(!ifs)
	{
		cout << "Cannot find the ILI file: " << INPUT.ili_file << endl;
		exit(0);
	}
	else
	{
		cout << "Open the ILI file: " << INPUT.ili_file << endl;
	}
	
	// output data
	ofstream ofs2("Joint_conditional_dist.dat");
	ofstream ofs3("Pos_ILI.dat");
	ofstream ofs4("ION.dat");
	ofstream ofs5("WATER.dat");


//-----------------  CORE CODE -------------------------------------

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "snapshot" << setw(12) << igeo << endl;

		proximity(ifs, cel, ofs3);
	}	

//------------------------------------------------------------------
	
	assert(count_geometry_number>0);
	double ave_na = 0.0;
	for(int ir=0; ir<nmesh; ++ir)
	{
		ave_na += record_water[0][ir];
	}
	ave_na/=(double)count_geometry_number;
	ofs_running << "average counted number of oxygen atoms = " << ave_na << endl;


	if(INPUT.only_hydroxide==true) ofs4 << "x O Oacc Odon H avgD avgA" << endl;
	ofs5 << "x O Odon0 Odon1 Odon2 Oacc0 Oacc1 Oacc2 avgD avgA avgTot" << endl;

	// 96 water molecules in a 12.445*12.445*18.5162 Angstroms^3 cell 
	// has a density of 1.0 g/cm3
	// 96/12.445/12.445/18.5162 = 0.0334756923 number_of_O/Angstrom^3
	const double fac = (double)count_geometry_number*(INPUT.celldm1*INPUT.celldm2*dr)*0.03347569;
	assert(fac!=0.0);
	ofs_running << "factor is " << fac << endl;

	for(int ir=0; ir<nmesh; ++ir)
	{
		double xx = ir*dr+INPUT.mdp0;

		// print out the ion information
		ofs4 << xx << " " << record_ion[0][ir]/fac << " " << record_ion[1][ir]/fac 
		<< " " << record_ion[2][ir]/fac << " " << record_ion[3][ir]/fac; 
		if(record_ion[6][ir]>0)
		{
			ofs4 << " " << record_ion[4][ir]/record_ion[6][ir]
				<< " " << record_ion[5][ir]/record_ion[6][ir] << endl;
		}
		else
		{
			ofs4 << " 0 0" << endl;
		}

		// print out the water molecules information
		ofs5 << xx << " " << record_water[0][ir]/fac << " " << record_water[1][ir]/fac
		<< " " << record_water[2][ir]/fac << " " << record_water[3][ir]/fac
		<< " " << record_water[4][ir]/fac << " " << record_water[5][ir]/fac
		<< " " << record_water[6][ir]/fac;

		if(record_water[9][ir]>0)
		{
			ofs5 << " " << record_water[7][ir]/record_water[9][ir] 
				<< " " << record_water[8][ir]/record_water[9][ir] 
				<< " " << (record_water[7][ir]+record_water[8][ir])/record_water[9][ir]
				<< endl;
		}
		else
		{
			ofs5 << " 0 0 0" << endl;
		}
	}


	// print out the joint conditional distribution function
	double sum=0.0;
	for(int ix=0; ix<INPUT.u1; ++ix)
	{
		for(int iy=0; iy<INPUT.u2; ++iy)
		{
			sum += joint_cd[ix][iy];
		}
	}
	// cout << "sum joint_cd is " << sum << endl;
	if(sum>0.0)
	{
		for(int ix=0; ix<INPUT.u1; ++ix)
		{
			for(int iy=0; iy<INPUT.u2; ++iy)
			{
				ofs2 << joint_cd[ix][iy]/sum << " ";	
			}
			ofs2 << endl;
		}
	}


	ofs_running << "proximity between " << INPUT.a0 << " " << INPUT.a1 << endl;
	ofs_running << "HBs of OH- is not included." << endl;
	if(ave_count>0)
	{
		ave_donate = ave_donate/ave_count;
		ave_accept = ave_accept/ave_count;
		ofs_running << "average donate HBs number " << ave_donate << endl;
		ofs_running << "average accept HBs number " << ave_accept << endl;
	}


	//-------------
	// clean up
	//-------------
	// close files
	ifs.close();
	ofs2.close();
	ofs3.close();
	ofs4.close();
	ofs5.close();

	for(int i=0; i<7; ++i) delete[] record_ion[i];
	delete[] record_ion;

	for(int i=0; i<10; ++i) delete[] record_water[i];
	delete[] record_water;

	// delete joint_cd
	for(int i=0; i<INPUT.u1; ++i)
	{
		delete[] joint_cd[i];
	}
	delete[] joint_cd;

	// delete arrays
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		delete[] interface[ix];
	}
	delete[] interface;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			delete[] gradient[ix][iy];
		}
		delete[] gradient[ix];
	}
	delete[] gradient;


	cout << "All done. Have a great day!" << endl;
	ofs_running << "All done. Have a great day!" << endl;

	return;
}


void MDP::proximity(ifstream &ifs, const Cell &cel, ofstream &ofs3)
{
	TITLE("MDP","proximity");

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

	// read the Instantaneous Liquid Interface (ILI)
	MDP::read_ili(ifs, this->interface, this->gradient);

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);

	// record the distance from each ion to the ILI interface
	double** aI_dis = new double*[INPUT.ntype];

	// print out the distance of each atom to the ILI
	ofs3 << cel.snapshot_index << " " << cel.snapshot_time << endl;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		aI_dis[it] = new double[cel.atom[it].na];
		for(int ia=0; ia<cel.atom[it].na;++ia)
		{
			double posx = cel.atom[it].pos[ia].x;
			double posy = cel.atom[it].pos[ia].y;
			double posz = cel.atom[it].pos[ia].z;
			while( posx >= norm1 ) posx -= norm1;	
			while( posx < 0 ) posx += norm1;
			while( posy >= norm2 ) posy -= norm2;	
			while( posy < 0 ) posy += norm2;
			if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide") 
			{
				while( posz >= norm3 ) posz -= norm3;	
				while( posz <0 ) posz += norm3;
			}
			else while( posz > INPUT.upper_z ) posz -= norm3;	
			int six=-1; int siy=-1;
			double drx=0.0; double dry=0.0; double drz=0.0;
			which_surface(six, siy, norm1, norm2, posx, posy, posz, drx, dry, drz, this->interface);
			double aaa = -(drx * gradient[six][siy][0] + dry * gradient[six][siy][1] + drz * gradient[six][siy][2]);

			// distance of all atoms to ILI
			ofs3 << cel.atom[it].id << " " << ia+1 << " " << aaa << endl;
			aI_dis[it][ia]=aaa;

			if(it == ito)
			{
				if(INPUT.nacc!=-1 and water[ia].naccept != INPUT.nacc) continue;
				double mdp0=INPUT.mdp0;

				if(aaa<mdp0) continue;
				int index = int((aaa-mdp0)/this->dr);
				if(index<0 or index>=this->nmesh) continue;
				
				record_water[0][index]+=1.0;

				if(water[ia].nH==2)
				{
					if(water[ia].ndonate==0)      record_water[1][index]+=1.0;
					else if(water[ia].ndonate==1) record_water[2][index]+=1.0;
					else if(water[ia].ndonate==2) record_water[3][index]+=1.0;
					if(water[ia].naccept==0)      record_water[4][index]+=1.0;
					else if(water[ia].naccept==1) record_water[5][index]+=1.0;
					else if(water[ia].naccept==2) record_water[6][index]+=1.0;

					record_water[7][index] += (double)water[ia].ndonate; 
					record_water[8][index] += (double)water[ia].naccept; 
					record_water[9][index] += 1.0; 
				}

				// calculate the joint conditional distribution funciton
				if(aaa>INPUT.a0 and aaa<INPUT.a1)
				{
					joint_conditional_distribution(aaa, ith, cel, posx, posy, posz, gradient[six][siy]);
					if(water[ia].nH==2)
					{
						ave_accept += water[ia].naccept;
						ave_donate += water[ia].ndonate;
						ave_count += 1;
					}
				}

			} // end it==ito
			
		}// end ia
	}// end it

	// second part to record the ion information
	if(INPUT.only_hydroxide==true)
	{
//		ofs_running << "water::nions " << Water::nions << endl;
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if(water[ia].nH==2) continue;

			int index=-1;
			double aaao=0.0;
			double aaah=0.0;
			int iao=-1;
			int iah=-1;

			aaao = aI_dis[ito][ia];
			ofs_running << "ion_position " << cel.snapshot_time << " " << aaao << " " << water[ia].naccept << endl;

			index = int((aaao-INPUT.mdp0)/dr);
			if(index<0 or index>=nmesh) continue;
			record_ion[0][index]++;

			for(int ia1=0; ia1<water[ia].naccept; ++ia1)
			{
				iao=water[ia].acceptO[ia1];
				aaao = aI_dis[ito][iao];
				index = int((aaao-INPUT.mdp0)/dr);
				if(index<0 or index>=nmesh) continue;
				record_ion[1][index]++;
			} 
			for(int ia2=0; ia2<water[ia].ndonate; ++ia2)
			{
				iao=water[ia].donateO[ia2];
				aaao = aI_dis[ito][iao];
				index = int((aaao-INPUT.mdp0)/dr);
				if(index<0 or index>=nmesh) continue;
				record_ion[2][index]++;
			}
			for(int ia3=0; ia3<water[ia].nH; ++ia3)
			{
				iah=water[ia].indexH[ia3];
				aaah = aI_dis[ith][iah];
				index = int((aaah-INPUT.mdp0)/dr);
				if(index<0 or index>=nmesh) continue;
				record_ion[3][index]++;
			}
			record_ion[4][index] += (double)water[ia].ndonate; 
			record_ion[5][index] += (double)water[ia].naccept; 
			record_ion[6][index] += 1.0; 

		}
	}

	// delete water and aI_dis;
	delete[] water;

	for(int it=0; it<INPUT.ntype; ++it)
	{
		delete[] aI_dis[it];
	}
	delete[] aI_dis;

	return;
}

void MDP::read_ili(ifstream &ifs, double** interface_in, double*** grad_in)
{
	string title;
	READ_VALUE(ifs, title);
	//cout << title << endl;
	int nx, ny;
	ifs >> nx >> ny;
	assert(INPUT.nx == nx);
	assert(INPUT.ny == ny);
	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			ifs >> interface_in[ix][iy];
		}
	}

	for(int j=0; j<3; ++j)
	{
		READ_VALUE(ifs,title);
		ifs >> nx >> ny;
		for(int ix=0; ix<nx; ++ix)
		{
			for(int iy=0; iy<ny; ++iy)
			{
				ifs >> grad_in[ix][iy][j];
			}
		}
	}

	// unit 
	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			double sum = sqrt(grad_in[ix][iy][0]*grad_in[ix][iy][0]+
			grad_in[ix][iy][1]*grad_in[ix][iy][1]+
			grad_in[ix][iy][2]*grad_in[ix][iy][2]);
	
			if(sum>0.0)
			{				
				grad_in[ix][iy][0]/=sum;
				grad_in[ix][iy][1]/=sum;
				grad_in[ix][iy][2]/=sum;
			}
		}
	}


	return;
}

void MDP::which_surface(int &six, int &siy,
const double &norm1, const double &norm2,
const double &posx, const double &posy, const double &posz, 
double& drx, double& dry, double& drz, double** interface_in)
{
	// calculate which surface point this water belongs to
	six=-1; siy=-1;
	drx = 0.0; dry = 0.0; drz = 0.0;
	double maxdis = 100000;
	const double dx = norm1/(double)INPUT.nx;
	const double dy = norm2/(double)INPUT.ny;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			double x = dx*ix;	
			double y = dy*iy;
			double dx = shortest(x, posx, norm1);
			double dy = shortest(y, posy, norm2);
			double dz = interface_in[ix][iy]-posz; 
			double dis = sqrt(dx*dx + dy*dy + dz*dz);
			if(dis < maxdis)
			{
				maxdis = dis;
				six = ix; siy = iy;
				drx = dx; dry = dy; drz = dz;
			}
		}
	}
//	cout << maxdis << " " << six << " " << siy << endl;
	assert(six!=-1); assert(siy!=-1);
	return;
}


// compute the joint conditional distribution function
void MDP::joint_conditional_distribution(const double &aaa, const int &ith, const Cell &cel, 
const double& posx, const double& posy, const double& posz, const double *grad)
{
	// compute H neighbours of oxygen
	int countH=0;
	int iu1=0;
	int iu2=0;
	for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
	{
		double ho_dx = shortest(posx, cel.atom[ith].pos[ia2].x, INPUT.celldm1);
		double ho_dy = shortest(posy, cel.atom[ith].pos[ia2].y, INPUT.celldm2);
		double ho_dz = shortest(posz, cel.atom[ith].pos[ia2].z, INPUT.celldm3);
		double ho_dis = sqrt(ho_dx*ho_dx+ho_dy*ho_dy+ho_dz*ho_dz);
		if(ho_dis<INPUT.rcut_oh)
		{
			++countH;
			if(countH==1)
			{
				double ucos1=(ho_dx*grad[0]+ho_dy*grad[1]+ho_dz*grad[2])/ho_dis;
				iu1 = (1+ucos1)/du1;
				//cout << "ucos1=" << ucos1 << " iu1=" << iu1 << endl;
			}
			else if(countH==2)
			{
				double ucos2=(ho_dx*grad[0]+ho_dy*grad[1]+ho_dz*grad[2])/ho_dis;
				iu2 = (1+ucos2)/du2;
				//cout << "ucos2=" << ucos2 << " iu2=" << iu2 << endl;
				joint_cd[iu1][iu2]+=1;
				break;
			}

			// mohan added 2018-08-12

		}
	}
	return;
}

