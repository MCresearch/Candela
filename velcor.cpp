#include "cellFile.h"
#include "input.h"
#include "velcor.h"

void VelCor::Routine()
{
	TITLE("VelCor","Routine");

	cal(); 

	return;
}


void VelCor::cal()
{
	TITLE("VelCor","cal");

	// --> INITIALIZE <--
	
	// (1.1) total number of atoms.
	const int nat = INPUT.natom;
	const int nat_used = INPUT.velcor_atom; 
	assert(INPUT.velcor_atom>0);
	cout << " Number of atoms used to calculate velocity correlation functions : " << nat_used << endl;

	// (1.2)* neqi = number of equilibrium states.
	// we need 'neqi' configurations to sample the ensemble.
	// INPUT.velcor_1=1 means the calculation starts from ion.1.dat
	const int neqi = INPUT.velcor_neqi; 
	cout << " Number of configurations used to represent the ensemble : " << neqi << endl;
	
	// length of correlation we want to calculate.
	// number of time steps, for example 1~100 = (100-1) + 1 = 100 steps.
	const int nt = INPUT.velcor_2-INPUT.velcor_1+1;
	cout << " Length of correlation we want to calculate, configureation number : " << nt << endl;

	assert(nat>0);
	assert(neqi>0);
	assert(nt>0);

	bool *file_exist = new bool[nt];
	bool *file_exist0 = new bool[neqi];

	// vc = velocity correlation function
	// set this to float to increase the speed.
	float** vc = new float*[nat];
	for(int ia=0; ia<nat; ++ia)
	{
		vc[ia] = new float[nt];
		ZEROS(vc[ia], nt);
	}

	// step interval for dynamics 
	const int step_interval_dynamics= INPUT.step_interval_dynamics;
	cout << " To prevent the correltaions, we use configurations every "
	 << step_interval_dynamics << " steps." << endl;
	assert(step_interval_dynamics>0);




	// --> BODY <--

	// how many atoms we want to calculate the 
	// velocity autocorrelation function.
	for(int ia=0; ia<nat_used; ++ia)
	{
		//ia=43;//mohan bad
		cout << " Atom " << ia+1 << endl;
		int count = 0;
		ZEROS(file_exist0, neqi);

		//// How many eqilibrium states we choose to represent sesemble,
		//// the more we choose, the more smooth the curves are.
		for(int iequi=0; iequi<neqi; ++iequi)
		{
			CellFile cel1;
			stringstream ss;

			// step interval for dynamics 
			// for exmaple, if the MD time interval is 0.1fs,
			// 250 * 0.1fs = 25fs,
			// we choose configuration between every 25 fs,
			// if we choose 200 configurations for ensemble,
			// we need 5 ps of simulation.
			const int file_index = iequi * step_interval_dynamics + INPUT.velcor_1;
			cout << " For ensemble " << iequi+1 <<", should read in file : ion." << file_index << ".dat" << endl;
			
			// try to find if the file exists. 
			// if file doesn't exist, go on for next.
			ss << file_index;
			cel1.file_name=ss.str();

			file_exist0[iequi]=CellFile::ReadVelocity(cel1);
			if(!file_exist0[iequi]) continue;

			// calculate <v1(0)*v1(0)>
			// velcor_atom-1 becauses the first atom starts from 0 index.
			const float v1v1 = this->velocity_correlation_functions(cel1,cel1,ia);
			// cout << " <v1 | v1>=" << v1v1 << endl;

			ZEROS(file_exist, nt);
			// The length of correlation functions are represented by the 
			// steps 'nt'. Typical value is 3ps, 30000 steps for MD time
			// interval 0.1fs, if each file represent 1fs, I need 3000
			// files.
			// cout << " nt=" << nt << endl;
			for(int iv=0; iv<nt; ++iv)
			{
				const int ifile = iv+file_index;
				// cel_in : input geometry file
				CellFile cel2;
				stringstream ss;
				ss << ifile;// Starts from 1.
				cel2.file_name=ss.str();
				//cout << " File name is " << cel2.file_name << endl;
				file_exist[iv] = CellFile::ReadVelocity(cel2);
				//// if the file doesn't exist, continue.
				if(!file_exist[iv]) continue;
				else
				{
				//	cout << cel2.file_name << " Exist." << endl;
				}
				//// calculate the correlation function
				//// <v1(t)*v1(0)>
				float this_vc=this->velocity_correlation_functions(cel1,cel2,ia);

				//block this for tmp.
				//this_vc/=v1v1;
				vc[ia][iv]+=this_vc;
				//cout << " vc[" << iv << "]=" << vc[ia][iv] << endl;
			}
			++count;
		}// end iequi

		cout << " Ensemble number for VAF for atom " << ia+1 << " : "  << count << endl;

		assert(count>0);
		for(int iv=0; iv<nt; ++iv)
		{
			vc[ia][iv]/=count;
		}
	}// end ia

	this->write_vc(vc, nt, file_exist);

	// --> CLEAN <--
	delete[] file_exist;
	delete[] file_exist0;

	for(int ia=0; ia<nat; ++ia)
	{
		delete[] vc[ia];
	}
	delete[] vc;

	return;
}

float VelCor::velocity_correlation_functions(
	const Cell &cel1,
	const Cell &cel2,
	const int &ia
)
{
	TITLE("VelCor","velocity_correlation_function");

	int it=0;
	float product = 0.0;

	// calculate ( v1 \cdot v2 )
	product = cel1.atom[it].vel[ia].x * cel2.atom[it].vel[ia].x
	+ cel1.atom[it].vel[ia].y * cel2.atom[it].vel[ia].y
	+ cel1.atom[it].vel[ia].z * cel2.atom[it].vel[ia].z;

//	cout << " cel1 vel=" << cel1.atom[it].vel[ia].x << " " << cel1.atom[it].vel[ia].y << " " << cel1.atom[it].vel[ia].z << endl;
//	cout << " cel2 vel=" << cel2.atom[it].vel[ia].x << " " << cel2.atom[it].vel[ia].y << " " << cel2.atom[it].vel[ia].z << endl;
//	cout << " Product =" << product << endl;

	return product;
}

void VelCor::write_vc(
	float** vc,
	const int &nt, 
	const bool* file_exist) const
{
	TITLE("VelCor","write_vc");

	ofstream ofs(INPUT.velcor_out.c_str());
	cout << " Output the velocity correlation function." << endl;
	int i=1;
	for(int iv=0; iv<nt; ++iv)
	{
		if(file_exist[iv])
		{
			float sum_vc = 0.0;
			for(int ia=0; ia<INPUT.natom; ++ia)
			{
				sum_vc += vc[ia][iv];
			}
			sum_vc /= (float)INPUT.velcor_atom;
			ofs << i << " " << sum_vc << endl;
			//ofs << i << " " << vc[0][iv] << " " << vc[1][iv] << endl;
			++i;
		}
	}
	// cout << " i=" << i << " nt=" << nt << endl;
	ofs.close();
	return;
}


