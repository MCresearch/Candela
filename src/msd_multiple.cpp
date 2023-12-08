#include "cellFile.h"
#include "input.h"
#include "msd_multiple.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"

MSD_Single::MSD_Single()
{
	count_msd = 0;
	natom = 0;
	sx = sy = sz = 0.0;
	mx = my = mz = 0.0;
	saved_ion = -1;
}
MSD_Single::~MSD_Single(){}

void MSD_Single::allocate(const double &t0_in)
{
	assert(INPUT.msd_dt>0.0);
	assert(INPUT.msd_t>0);
	this->ndim = INPUT.msd_t / INPUT.msd_dt;
	this->msd = new double[ndim]();
	this->mmm = new int[ndim]();
	this->t = new double[ndim]();

	this->t0 = t0_in;
	this->t1 = t0_in + INPUT.msd_t;

	ofs_running << "t0 is " << t0 << " t1 is " << t1 << endl;
	for(int it=0; it<ndim; ++it)
	{
		t[it] = it*INPUT.msd_dt+t0;
		msd[it] = 0; // renxi added 20200505
	}
	assert(INPUT.msd_natom>0);
	this->pre_wpos = new Vector3<double>[INPUT.msd_natom];
	this->wpos = new Vector3<double>[INPUT.msd_natom];
	this->wpos0 = new Vector3<double>[INPUT.msd_natom];

	return;
}

void MSD_Single::deallocate()
{
	delete[] msd;
	delete[] mmm;
	delete[] t;
	delete[] pre_wpos;
	delete[] wpos;
	delete[] wpos0;
}

MSD_Multiple::MSD_Multiple()
{
}

MSD_Multiple::~MSD_Multiple(){}

void MSD_Multiple::Routine()
{
	TITLE("MSD_Multiple","Routine");
	
	if(RANK==0) cout << " Compute Multiple Mean Square Displacements" << endl;
	
	// the input parameters needed by computing multiple MSD. 
	ofs_running << " msd_n  = " << INPUT.msd_n << endl;
	ofs_running << " msd_t0 = " << INPUT.msd_t0 << endl;
	ofs_running << " msd_t  = " << INPUT.msd_t << endl;
	ofs_running << " msd_dt0 = " << INPUT.msd_dt0 << endl;
	ofs_running << " msd_dt = " << INPUT.msd_dt << endl;
	ofs_running << " msd_natom = " << INPUT.msd_natom << endl;
	ofs_running << " msd_stokes = " << INPUT.msd_stokes << endl;

	assert(INPUT.msd_n > 0);
	assert(INPUT.msd_t0 >= 0.0);
	assert(INPUT.msd_t > 0.0);
	assert(INPUT.msd_dt > 0.0);
	assert(INPUT.msd_natom > 0);

	// Each MSD is a member of MSD_Single
	ms = new MSD_Single[ INPUT.msd_n ];	
	for(int i=0; i<INPUT.msd_n; ++i)
	{
		// the input is the starting point
		ms[i].allocate(INPUT.msd_t0+INPUT.msd_dt0*i);
	}

	// open the output file
	ofstream ofs_msd_each;
	ofstream ofs_msd_total;

	if(RANK==0)
	{
		ofs_msd_each.open("MSD_each.txt");
    	ofs_msd_total.open("MSD_total.txt");
	}
	
	//**************************
	// BEGIN CALCULATING DATA
	//**************************
	int ipt=0;
	int geo_count=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
		//cout << "read geometry end." << endl;
		if(cel.read_and_used==false) 
		{
			cel.clean(); // mohan added 2018-05-27
			continue;
		}
		if(RANK==0) cout << "igeo=" << igeo << " " << cel.snapshot_time << endl;
		ofs_running << "igeo=" << igeo << endl;
		compute_msd(cel, igeo);
		geo_count++;
		cout << "geo_count = " << geo_count << endl;
		cel.clean(); // mohan added 2018-05-27, very important
#ifdef __MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}	

	cout << "total number of geometry file: " << geo_count << endl;

#ifdef __MPI
	for(int i=0; i<INPUT.msd_n; ++i)
	{
		// reduce 'msd'
		int dim=ms[i].ndim;
		double* tmp = new double[dim];
		for(int j=0; j<dim; ++j)
		{
			tmp[j] = ms[i].msd[j];
		}
	    MPI_Allreduce(tmp, ms[i].msd, dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		delete[] tmp;
		
		// reduce 'mmm'
		int* tmp2 = new int[dim];
		for(int j=0; j<dim; ++j)
		{
			tmp2[j] = ms[i].mmm[j];
		}
		MPI_Allreduce(tmp2, ms[i].mmm, dim, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		delete[] tmp2; 
	}
#endif


	// get the average msd	
	for(int i=0; i<ms[0].ndim; ++i)
	{
		for(int j=0; j<INPUT.msd_n; ++j)
		{
			if(ms[j].mmm[i]>0)
			{
				 ms[j].msd[i]/=(double)ms[j].mmm[i];
			}
		}
	}

	//*******************
	// print out the data
	//*******************
	if(RANK==0)
	{
		for(int i=0; i<ms[0].ndim; ++i)
		{
			// only need to print out the time of first set of data 
			double ttt = ms[0].t[i] - INPUT.msd_t0;
			if(ttt==0.0) continue;

			ofs_msd_each << ttt;
			ofs_msd_total << ttt;

			for(int j=0; j<INPUT.msd_n; ++j)
			{
				ofs_msd_each << " " << ms[j].msd[i];
			}

			double average=0.0;
			int nnn=0;
			for(int j=0; j<INPUT.msd_n; ++j)
			{
				if(ms[j].msd[i]>0.0)
				{
					average += ms[j].msd[i];
					++nnn;
				}
			}

			assert(nnn>=0);		
			if(nnn>0)
			{
				average/=nnn;
				ofs_msd_total << " " << average << " " << average/6.0/ttt;
			}
			//else if(nnn==0)
			//{
			//	ofs_msd << " 0";
			//}

			// compute the standard deviation
			if(INPUT.msd_n>0)
			{
				double sd=0.0;
				for(int j=0; j<INPUT.msd_n; ++j)
				{
					if(ms[j].msd[i]>0.0)
					{
						sd += (ms[j].msd[i]-average)*(ms[j].msd[i]-average);
					}
				}
				if(nnn>0)
				{
					if(nnn>1)
					{
						sd /= (nnn-1);
						sd = sqrt(sd);
						ofs_msd_total << " " << sd;
					}
					else
					{
						ofs_msd_total << " 0";
					}
				}
			}	

			ofs_msd_total << endl;
			ofs_msd_each << endl;
		}

		ofs_msd_each.close();
		ofs_msd_total.close();
	}

	// clean up
	for(int i=0; i<INPUT.msd_n; ++i)
	{
		ms[i].deallocate();
	}

	delete[] ms;

	return;
}


void MSD_Multiple::compute_msd(const Cell &cel, const int &igeo)
{
	// get ito, ith, and itc.
	int ito=-1;
	int ith=-1;
	int itc=-1;
	int it_select = -1;

	int na=0;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
		na += cel.atom[it].na;

		if (INPUT.ele_select != "none" and cel.atom[it].id == INPUT.ele_select)
		{
			it_select = it;
		}
	}
	//if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
	//if(INPUT.ntype==3){ assert(itc>=0); }

	if (INPUT.system == "water" or INPUT.system == "hydroxide" or INPUT.system == "hydronium")
	{
		assert(ito>=0);
		assert(ith>=0);
	}

	assert(INPUT.celldm1>0);
	assert(INPUT.celldm2>0);
	assert(INPUT.celldm3>0);

	Water *water;

	if(INPUT.system=="water" || INPUT.system=="hydronium" || INPUT.system=="hydroxide")
	{ 
		water = new Water[cel.atom[ito].na];
		Water::nions=0;
		HBs::setup_water(cel, water);
	}

	for(int i=0; i<INPUT.msd_n; ++i)
	{
//		ofs_running << "i= " << i << endl;
		each_msd(cel, ito, it_select, ms[i], water);
	}

	if(INPUT.system=="water" || INPUT.system=="hydronium" || INPUT.system=="hydroxide")
	{	
		delete[] water;
	}
	return;
}

void MSD_Multiple::each_msd(const Cell &cel, const int &ito, const int &it_select, MSD_Single &ms, const Water* water)
{
	if(cel.snapshot_time < ms.t0) return;
	else if(cel.snapshot_time>= ms.t1) return;

    // step 1
	bool flag_Stokes=INPUT.msd_stokes;
	int ion_index=-1;

	// the number of ions in solution should be 1 without any ambiguity
	if(INPUT.system=="hydronium" or INPUT.system=="hydroxide")
	{
		if(Water::nions!=1) return;
	}

	// record the jumps due to proton transfer events
	// works only for ions
	if(flag_Stokes==true)
	{

		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			if(INPUT.system=="hydronium" and water[ia].nH!=3) continue;
			else if(INPUT.system=="hydroxide" and water[ia].nH!=1) continue;
			ion_index=ia; // get the ion index
			break;
		}

		//ofs_running << "saved_ion=" << ms.saved_ion << endl;

		if(ms.saved_ion==-1)
		{
			ms.saved_ion=ion_index;
			ms.mx = ms.my = ms.mz = 0.0;
		}
		else if(ion_index!=ms.saved_ion)
		{
			ms.saved_ion=ion_index;
			ofs_running << "protontransfer_jumps " << cel.snapshot_index << " " << cel.snapshot_time << " " << ion_index+1;
			ms.mx += shortest(ms.sx, cel.atom[ito].pos[ion_index].x, INPUT.celldm1);
			ms.my += shortest(ms.sy, cel.atom[ito].pos[ion_index].y, INPUT.celldm2);
			ms.mz += shortest(ms.sz, cel.atom[ito].pos[ion_index].z, INPUT.celldm3);
			ofs_running << " mxyz: " << ms.mx << " " << ms.my << " " << ms.mz << endl;
		}

		// the position of the saved ion
		ms.sx = cel.atom[ito].pos[ion_index].x;
		ms.sy = cel.atom[ito].pos[ion_index].y;
		ms.sz = cel.atom[ito].pos[ion_index].z;
	}

	// step 1.2 compute the mass center//qianrui delete this part
/*
	double xmc=0.0;	
	double ymc=0.0;	
	double zmc=0.0;	
	double sum_mass=0.0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			double xx=cel.atom[it].pos[ia].x;
			double yy=cel.atom[it].pos[ia].y;
			double zz=cel.atom[it].pos[ia].z;
			while(xx<0) xx+=INPUT.celldm1;	
			while(xx>=INPUT.celldm1) xx-=INPUT.celldm1;	
			while(yy<0) yy+=INPUT.celldm2;	
			while(yy>=INPUT.celldm2) yy-=INPUT.celldm2;	
			while(zz<0) zz+=INPUT.celldm3;	
			while(zz>=INPUT.celldm3) zz-=INPUT.celldm3;	
			xmc += xx*cel.atom[it].mass;
			ymc += yy*cel.atom[it].mass;
			zmc += zz*cel.atom[it].mass;
			sum_mass += cel.atom[it].mass;
		}
	}
	assert(sum_mass>0.0);
	xmc /= sum_mass;
	ymc /= sum_mass;
	zmc /= sum_mass;	
//	ofs_running << "mass_center " << xmc << " " << ymc << " " << zmc << endl;
	
*/
	// step 2
	if(ms.count_msd==0)
	{
		int iat=0;
		int iaq=0;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			if (it_select != -1 and it != it_select)
			{
				continue;
			}
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				if(INPUT.system=="water")
				{
					//nothing happens	
					iaq = iat;
				}
				else if(INPUT.system=="hydronium" || INPUT.system=="hydroxide") 
				{
					if(it!=ito) continue;
					if(INPUT.system=="hydronium" and water[ia].nH!=3) continue;
					else if(INPUT.system=="hydroxide" and water[ia].nH!=1) continue;
				}
				else
				{
					iaq = iat;
				}
				ms.wpos[iaq] = cel.atom[it].pos[ia];
				ms.pre_wpos[iaq] = cel.atom[it].pos[ia];
				ms.wpos0[iaq] = cel.atom[it].pos[ia]; 
				/*qianrui delete it
				ms.xmc0=xmc;
				ms.ymc0=ymc;
				ms.zmc0=zmc;*/
				++iat;
			}
		}
	}
	else
	{
		Vector3<double> move;
		// mohan fixed a bug 2018-05-28, 0.145/0.005=28 is wrong
		double error = 0.0001; // (ps)
		//int index = (cel.snapshot_time - ms.t0 + error) / INPUT.msd_dt;
		int index = round((cel.snapshot_time - ms.t0) / INPUT.msd_dt);
		//cout << "time=" << cel.snapshot_time << " dt=" << INPUT.msd_dt << " t0=" << ms.t0 << endl; // test
		assert(index>=0);
		if (index >= ms.ndim) return;
		assert(index<ms.ndim);
		int iat=0;
		int iaq=0;
//		cout << "msd0=" << ms.msd[index] << " mmm0=" << ms.mmm[index] << " index=" << index << endl;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			if (it_select != -1 and it != it_select)
			{
				continue;
			}
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
#ifdef __MPI
				if( iat%NPROC!=RANK ) {iat++; continue;}
//				ofs_running << setw(10) << "Atom" << setw(10) << iat << endl;
#endif

				if(INPUT.system=="water")
				{
					//nothing happens	
					iaq = iat;
				}
				else if(INPUT.system=="hydronium" || INPUT.system=="hydroxide") 
				{
					if(it!=ito) continue;
					if(INPUT.system=="hydronium" and water[ia].nH!=3) continue;
					else if(INPUT.system=="hydroxide" and water[ia].nH!=1) continue;
				}
				else
				{
					iaq = iat;
				}
				
				move.x = shortest(ms.pre_wpos[iaq].x, cel.atom[it].pos[ia].x, INPUT.celldm1);
				move.y = shortest(ms.pre_wpos[iaq].y, cel.atom[it].pos[ia].y, INPUT.celldm2);
				move.z = shortest(ms.pre_wpos[iaq].z, cel.atom[it].pos[ia].z, INPUT.celldm3);
				ms.wpos[iaq] = ms.wpos[iaq] - move;
				//ms.pre_wpos[iaq] = cel.atom[it].pos[ia];
				double dx = ms.wpos[iaq].x-ms.wpos0[iaq].x;
				double dy = ms.wpos[iaq].y-ms.wpos0[iaq].y;
				double dz = ms.wpos[iaq].z-ms.wpos0[iaq].z;
				/*
				if( dx*dx + dy*dy + dz*dz == 0 )
				{
					cout << iaq << " " << dx*dx + dy*dy + dz*dz << endl; // test
					cout << ms.wpos0[iaq].x << " " << ms.wpos0[iaq].y << " " << ms.wpos0[iaq].z << endl;
					cout << dx << " " << dy << " " << dz << endl;
					cout << cel.atom[it].pos[ia].x << " " << cel.atom[it].pos[ia].y << " " << cel.atom[it].pos[ia].z << endl;
					cout << move.x << " " << move.y << " " << move.z << endl;
					cout << ms.pre_wpos[iaq].x << " " << ms.pre_wpos[iaq].y << " " << ms.pre_wpos[iaq].z << endl;
					cout << ms.wpos[iaq].x << " " <<  ms.wpos[iaq].y << " " <<  ms.wpos[iaq].z << endl;
					exit(0);
				}
				*/
				// renxi added 20200505
				ms.pre_wpos[iaq] = cel.atom[it].pos[ia];
				// mx, my, mz are due to the proton transfer jumps
				if(flag_Stokes==true)
				{
					dx += ms.mx;	
					dy += ms.my;	
					dz += ms.mz;	
				}

				// corrections due the movement of mass center
				/*qianrui delete it
				dx -= xmc - ms.xmc0;
				dy -= ymc - ms.ymc0;
				dz -= zmc - ms.zmc0;*/


				ms.msd[index] += dx*dx + dy*dy + dz*dz;  // unit is A^2
				ms.mmm[index] += 1;
				++iat;
				/*
				if(ms.msd[index]==0)
				{
					cout << dx << " " << dy << " " << dz << endl;
					exit(0);
				}
				*/
				// renxi added 20200507
			}
		}
//		cout << "msd=" << ms.msd[index] << " mmm=" << ms.mmm[index] << endl;
		if(INPUT.system=="hydronium" or INPUT.system=="hydroxide") assert(iat==1);
	}
	++ms.count_msd;

	return;
}

int MSD_Multiple::round(double r)
{
    if(r-floor(r)>=0.5){return ceil(r);}
    else 			   {return floor(r);}
}
