#include "wannier.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Wannier::Wannier() 
{
}

Wannier::~Wannier() 
{
}

void Wannier::infrared()
{
	this->tcor = INPUT.tcor;
	assert(tcor>0);
	this->dt_au=2.4188843E-05; // unit is ps
	// mohan update on 2018-07-05
	this->interval = INPUT.dt_snapshots/dt_au; // how many dt_au are separated by two snapshots
	
	ofs_running << " dt_snapshots = " << INPUT.dt_snapshots << endl;
	ofs_running << " interval = " << interval << endl;


	assert(interval>0);
	double *uacf = new double[tcor+1]();

	// generate the correlation function in "vUACF.dat"
	// and then do the Fourier transform
	if(INPUT.func_b==1)
	{
		correlation_function_vdipole(uacf);
		fourier_transform(uacf, tcor+1);
	}
	// read in the correlation function
	// and then do the Fourier transform
	else if(INPUT.func_b==2)
	{
		// read in the correlation function generated by func_b==1
		read_correlation_vdipole(uacf);

		// multiply the correlation function with a gaussian function in order to make the Fourier transform smooth
		double sigma=INPUT.factor;
		// cutoff for correlation function (need to observe the correlation function to determine the appropriate rcut)
		double tcut=INPUT.bdf_rcut; 
		ofstream ofs("smooth_vdipole.dat");
		cout << "tcut=" << tcut << endl;
		for(int it=0; it<=tcor; ++it)
		{
			double tt = dt_au*interval*(double)it;
			double aaa=-pow(tt-tcut,2.0)/(2*sigma*sigma);
			double bbb=exp(aaa);
			ofs << tt << " " << bbb << " " << uacf[it];
			if(tt>=tcut)
			{
				uacf[it] *= bbb; 
			}
			ofs << " " << uacf[it] << endl;
		}
		ofs.close();

		// perfrom the Fourier transfroms
		fourier_transform(uacf, tcor+1);
	}
/*
	else if(INPUT.func_b==3) // compute the dipole correlation function
	{
		correlation_function_dipole(uacf);
		fourier_transform(uacf, tcor+1);
	}
	else if(INPUT.func_b==4) // Fourier transform of dipole autocorrelation function
	{
		read_correlation_dipole(uacf);

		// multiply the correlation function with a gaussian function in order to make the Fourier transform smooth
		double sigma=INPUT.factor;
		// cutoff for correlation function (need to observe the correlation function to determine the appropriate rcut)
		double tcut=INPUT.bdf_rcut; 
		ofstream ofs("smooth_dipole.dat");
		cout << "tcut=" << tcut << endl;
		for(int it=0; it<=tcor; ++it)
		{
			double tt = dt_au*interval*(double)it;
			double aaa=-pow(tt-tcut,2.0)/(2*sigma*sigma);
			double bbb=exp(aaa);
			ofs << tt << " " << bbb << " " << uacf[it];
			if(tt>=tcut)
			{
				uacf[it] *= bbb; 
			}
			ofs << " " << uacf[it] << endl;
		}
		ofs.close();

		// perfrom the Fourier transfroms
		fourier_transform(uacf, tcor+1);
	}
*/

	delete[] uacf;
	return;
}

void Wannier::correlation_function_dipole(double* uacf)
{
	cout << "Compute infrared spectra based on dipole moments" << endl;
	
	int tend=INPUT.geo_2;
	int natom=INPUT.natom1; // number of molecules (should be oxygen in water)
	assert(natom>0);

	// open the file about dipoles
	ofstream out("UACF.dat");

	// file pointers
	FILE *inp1, *inp2;
	inp1 = fopen("all_dipole.dat","r");
	inp2 = fopen("all_dipole.dat","r");
	fpos_t position;

	double* norm1 = new double[tcor+1];
    for(int it=0;it<tcor+1;++it)
    {
        uacf[it]=0.0;
        norm1[it]=0.0;
    }

	// read in vdipole and calculate the correlation function
	double* ux0 = new double[natom]();
    double* uy0 = new double[natom]();
    double* uz0 = new double[natom]();
    double* uxt = new double[natom]();
    double* uyt = new double[natom]();
    double* uzt = new double[natom]();

	int snapshot=0;
	double time_now=0.0;
	double tdip=0.0;

	for(int t0=0; t0<tend; ++t0)
	{
		double tot_dx0=0.0;
		double tot_dy0=0.0;
		double tot_dz0=0.0;

		fscanf(inp1,"%d%lf",&snapshot,&time_now);
		cout << "snapshot " << snapshot << " time " << time_now << endl;

		for(int k=0;k<natom;++k)
		{
			fscanf(inp1,"%lf%lf%lf%lf",&ux0[k],&uy0[k],&uz0[k],&tdip);
			tot_dx0 += ux0[k];
			tot_dy0 += uy0[k];
			tot_dz0 += uz0[k];
		}

		uacf[0] += tot_dx0*tot_dx0+tot_dy0*tot_dy0+tot_dz0*tot_dz0; 
		norm1[0]++;

		fgetpos(inp1,&position);
		fsetpos(inp2,&position);

		int tt0max=0;
        if(tend<(t0+tcor))
        {
            tt0max=tend-t0;
        }
        else
        {
            tt0max=tcor;
        }
        printf("snapshot %d correlation with the following %d snapshots\n",t0, tt0max);

		if(tt0max==0) break;

		// for the following n snapshots
		for(int it=1;it<=tt0max;it++)
		{
			double tot_dxt=0.0;
			double tot_dyt=0.0;
			double tot_dzt=0.0;
			fscanf(inp2,"%d%lf",&snapshot,&time_now);
			for(int k=0;k<natom;++k)
			{
				fscanf(inp2,"%lf%lf%lf%lf\n",&uxt[k],&uyt[k],&uzt[k],&tdip);
				tot_dxt+=uxt[k];	
				tot_dyt+=uyt[k];	
				tot_dzt+=uzt[k];	
			}
            uacf[it] += tot_dx0*tot_dxt+tot_dy0*tot_dyt+tot_dz0*tot_dzt; 
            norm1[it]+= 1;
		}
	}
	
	// normalize, then output the final correlation function
    for(int it=0;it<tcor+1;++it)
    {
        if(norm1[it]!=0.0)
        {
           uacf[it]=uacf[it]/norm1[it];
           double tt=it*dt_au*interval; // time unit is ps
            out << tt << " " << uacf[it] << endl;
        }
        else
        {
            cout << "WARNING! norm1=0.0" << endl;
        }
    }

    delete[] ux0;
    delete[] uy0;
    delete[] uz0;
    delete[] uxt;
    delete[] uyt;
    delete[] uzt;
    delete[] norm1;

    // close
    fclose(inp1);
    fclose(inp2);

	return;
}

void Wannier::correlation_function_vdipole(double *uacf)
{
	cout << "Compute infrared spectra based on velocities of dipole moments" << endl;

	int tend=INPUT.geo_2;
	int natom=INPUT.natom1; // number of molecules (should be oxygen in water)
	assert(natom>0);

	// open the file about velocity of dipoles
	ifstream in(INPUT.vdipole_file.c_str());
	ifstream in1(INPUT.vdipole_file.c_str());	
	ifstream *in2;	
	ofstream out("vUACF.dat");

	// file pointers
	FILE *inp1, *inp2;
	inp1 = fopen(INPUT.vdipole_file.c_str(),"r");
	inp2 = fopen(INPUT.vdipole_file.c_str(),"r");
	fpos_t position;
	
	// calculate averaged dipole moments
    double totx=0.0;
    double toty=0.0;
    double totz=0.0;

   	for(int i=0; i<tend; ++i) 
    {
		for(int j=0; j<natom; ++j)
		{
			double dx,dy,dz;
			in >> dx >> dy >> dz;
        	totx += dx;
        	toty += dy;
        	totz += dz;
		}
    }

    double avgx = totx/tend/natom;
    double avgy = toty/tend/natom;
    double avgz = totz/tend/natom;


	ofs_running << "avgx=" << avgx << " avgy=" << avgy << " avgz=" << avgz << endl;


	double* norm1 = new double[tcor+1];
    for(int it=0;it<tcor+1;++it)
    {
        uacf[it]=0.0;
        norm1[it]=0.0;
    }

	// read in vdipole and calculate the correlation function
	double* ux0 = new double[natom]();
	double* uy0 = new double[natom]();
	double* uz0 = new double[natom]();
	double* uxt = new double[natom]();
	double* uyt = new double[natom]();
	double* uzt = new double[natom]();


	for(int t0=0; t0<tend; ++t0)
    {
        for(int k=0;k<natom;++k)
        {
		//	in1 >> ux0[k] >> uy0[k] >> uz0[k];
			fscanf(inp1,"%lf%lf%lf",&ux0[k],&uy0[k],&uz0[k]);
        }


		// compute the deviations of dipole velocities in terms of the average dipole velocity
		// at time 0, which is \sum_ij[ui(0)-uavg]*[uj(0)-uavg]
        for(int j=0;j<natom;++j)
        {
            for(int k=0;k<natom;++k)
            {
                uacf[0] += (ux0[j]-avgx)*(ux0[k]-avgx)+(uy0[j]-avgy)*(uy0[k]-avgy)+(uz0[j]-avgz)*(uz0[k]-avgz);
            }
        }
        norm1[0]+=1.0;

		// set the file pointer inp2 equals to the file pointer inp1 
		// in2 = &in1;
		fgetpos(inp1,&position);
		fsetpos(inp2,&position);

		// calculate the number of snapshots that can be used to compute the correlation function between dipoles 
		// (or velocities of dipoles)
		int tt0max=0;
        if(tend<(t0+tcor))
		{
            tt0max=tend-t0;
		}
        else 
		{
			tt0max=tcor;
		}
        //printf("snapshot %d correlation with the following %d snapshots\n",t0, tt0max);
		cout << "snapshot " << t0 << " is correlated with the following " << tt0max << "snapshots" << endl;
		

		if(tt0max==0) break;


		// for the following n snapshots 
        for(int it=1;it<=tt0max;it++)
        {
            for(int k=0;k<natom;++k)
            {
	//			*in2 >> uxt[k] >> uyt[k] >> uzt[k];
				fscanf(inp2,"%lf%lf%lf\n",&uxt[k],&uyt[k],&uzt[k]);
            }
            for(int j=0;j<natom;++j)
            {
				double xx = ux0[j]-avgx;
				double yy = uy0[j]-avgy;
				double zz = uz0[j]-avgz;
                for(int k=0;k<natom;++k)
                {
                    uacf[it] += xx*(uxt[k]-avgx)+yy*(uyt[k]-avgy)+zz*(uzt[k]-avgz);
                }
            }
			norm1[it]+= 1.0;
        }
    } // end reading


	// normalize, then output the final correlation function
    for(int it=0;it<tcor+1;++it)
    {
		if(norm1[it]!=0.0)
		{
     	   uacf[it]=uacf[it]/norm1[it];
     	   double tt=it*dt_au*interval; // time unit is ps 
			out << tt << " " << uacf[it] << endl;
		}
		else
		{
			cout << "WARNING!" << endl;
		}
    }

	cout << "avgx=" << avgx << " avgy=" << avgy << " avgz=" << avgz << endl;
	ofs_running << "avgx=" << avgx << " avgy=" << avgy << " avgz=" << avgz << endl;

	// close all files
	in.close();
	in1.close();

	delete[] ux0;
	delete[] uy0;
	delete[] uz0;
	delete[] uxt;
	delete[] uyt;
	delete[] uzt;
	delete[] norm1;

	// close 
	fclose(inp1);
	fclose(inp2);

	return;
}


void Wannier::read_correlation_dipole(double* uacf)
{
	ifstream ifs("UACF.dat");
	if(!ifs)
	{
		cout << "could not open file: UACF.dat" << endl;
	}
	else
	{
		cout << "reading the autocorrelation function." << endl;
	}

	double* tt = new double[tcor]();

	for(int i=0; i<=tcor; ++i)
	{
		ifs >> tt[i] >> uacf[i];
	}	

	delete[] tt;	

	ifs.close();
	cout << "Finish reading the autocorrelation function" << endl;


	ifs.close();
	return;
}


void Wannier::read_correlation_vdipole(double* uacf)
{
	ifstream ifs("vUACF.dat");
	if(!ifs)
	{
		cout << "could not open file: vUACF.dat" << endl;
		exit(0);
	}
	else
	{
		cout << "reading the autocorrelation function." << endl;
	}

	double* tt = new double[tcor]();

	for(int i=0; i<=tcor; ++i)
	{
		ifs >> tt[i] >> uacf[i];
	}	

	delete[] tt;	

	ifs.close();
	cout << "Finish reading the autocorrelation function" << endl;
		

	return;
}

void Wannier::fourier_transform(double* uacf, const int &ndim)
{
	assert(ndim>0);

	cout << "Perform Fourier transform" << endl;

	// open file: final infrared spectra
	ofstream ofs_in("infrared.dat");
	double* infrared = new double[ndim]();

	// fourier transfrom
	// time unit is ps
	double tmax = ndim*dt_au*interval; // unit is ps
	//double dw = 2*pi/tmax; // this is wrong!
	double dw = 1.0/tmax; // unit is 1/ps
	
	assert(INPUT.rho_ref>0);
	// unit is Angstrom^3
	assert(INPUT.natom1>0);
	double volume = INPUT.natom1 * 18 * 1.6605/ INPUT.rho_ref;
	volume *= 1.0e-30;
	cout << "volume is " << volume << " m^3" << endl;

	// temperature unit is K
	double temperature=INPUT.temperature;

	double beta = 1.0/(1.38064852*1.0e-23)/temperature;
	double c=299792458; // m/s

	// 1 Debye = 0.20819434 e*Angstrom
	// 1 e = 1.602*1.0e-19 C
	// change unit to C*m for M(0)
	double unit_basic = 0.20819434*1.602176565*1.0e-19*1.0e-10;
	// change unit to ps for dM(0)/dt
	double unit = unit_basic/2.418884326505*1.0e5;
	// because dot(M(0))*dot(M(t)) change unit to C^2 * m^2 / ps^2
	double unit2 = unit * unit;
	double unit_all = 2.0*PI*beta/3.0/c/volume*unit2;
	cout << "unit_all is " << unit_all << endl;




	// in the integral there is a dt, which is in ps,
	// but there are two t units in M(0) and M(t)
	double dt = 1.0e-12; 
	unit_all /= dt;


	// 1F = A^2 * s^4 * kg^-1 * m^-2
	double uc = unit_all/(4*PI*EPSILON0); // from Gs to SI unit, due to the C^2
	uc /= 100; // change m^-1 to cm^-1
    uc /= 1000; // change cm^-1 to 1000cm^-1

	cout << "uc=" << uc << endl;
	cout << "interval =" << interval << endl;
	cout << "dt_au = " << dt_au << endl;

	for(int iw=0; iw<ndim; ++iw)
	{
		infrared[iw]=0.0;
		double ww = dw*(double)iw;
		for(int it=0; it<ndim; ++it)
		{
			double tt = dt_au*interval*(double)it; // time unit is in ps
			//infrared[iw] += 2*uacf[it]*cos(ww*tt)*(dt_au*interval); // the dt_au*interval seems to be included in uacf already
			//infrared[iw] += 2*uacf[it]*cos(ww*tt);
			infrared[iw] += 2*uacf[it]*cos(ww*tt)*(dt_au*interval);
		}
		// in spectroscopy, "wavenumber" often refers to a frequency that has been
		// divided by the speed of light in vacuum, 0.01 is from 1/m to 1/cm-1
		ofs_in << ww*(1/(2*PI*c*1.0e-12)*0.01) << " " << infrared[iw]*uc << endl; 
	}
	delete[] infrared;
	ofs_in.close();

	return;
}

void Wannier::Routine()
{
	cout << "===== Compute dipole moments based on the MLWF Centers =====" << endl;	

	// print out all dipole moments for all water molecules 
	ofstream ofs_dipole("all_dipole.dat");
	// print out the dipole velocity
	ofstream ofs_vdipole("all_vdipole.dat");
	// print out the distribution of total dipoles
	ofstream ofs_tdipole("total_dipoleM.dat");
	ofs_dipole  << setprecision(8);
	ofs_vdipole << setprecision(8);
	ofs_tdipole << setprecision(8);
 	init_previous=false;	
	previous_t=-1.0;
	previous2_t=-1.0;

	// plot out the final distance between H and wannier centers
	this->dr=INPUT.dr;
	assert(INPUT.dr>0);
	double rcut=INPUT.rcut;
	this->nr=rcut/dr+10; //+10 for safety
	dis_Hwan = new double[nr]();

	// prepare for 2D plot
	// an example of input
	//  nx 70
	//  ny 70
	//  x0 0.95
	//  y0 0.25
	//  dx 0.0025
	//  dy 0.005
	/*
    assert(INPUT.nx>0);
    assert(INPUT.ny>0);
    int nx = INPUT.nx;
    int ny = INPUT.ny;
    this->coord_xy = new double*[nx];
    for(int ix=0; ix<nx; ++ix)
    {
        this->coord_xy[ix] = new double[ny]();
    }
	*/

	// prepare for dipole moments;
	assert(INPUT.dz>0.0);
	this->nd=INPUT.rcut1/INPUT.dz+10;
	this->dipole_dis = new double[nd]();

	this->avg_dipole = 0.0;	
	this->count_dipole=0;

	// for total dipole
	this->rcut_total_dipole = 100;
	this->dr_total_dipole = 0.2;
	int nmesh_total_dipole = rcut_total_dipole/dr_total_dipole+1;
	this->dis_total_dipole = new double[nmesh_total_dipole]();


	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

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


		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << cel.snapshot_index << " igeo " << igeo << endl;

		read_wan(cel, dis_Hwan, ofs_dipole, ofs_vdipole, ofs_tdipole);
		cel.clean();
	}	



	for(int i=0; i<nmesh_total_dipole; ++i)
	{
		ofs_tdipole << i * dr_total_dipole << " " << dis_total_dipole[i] << endl;
	}



	// print out dipole moments
	ofstream ofsd("distribution_dipole.dat");
	double sum=0.0;
	for(int id=0; id<nd; ++id)
	{
		sum += dipole_dis[id]*INPUT.dz;
	}
	for(int id=0; id<nd; ++id)
	{
		ofsd << id*INPUT.dz << " " << dipole_dis[id]/sum << endl;
	}	
	ofsd.close();
	delete[] dipole_dis;

	if(count_dipole>0)
	{
		ofs_running << "avg_dipole=" << avg_dipole/(double)count_dipole << endl;
		cout << "avg_dipole=" << avg_dipole/(double)count_dipole << endl;
	}

	// print out 1D Wannier centers
	ofstream ofs("distribution_MLWF.dat");

	sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_Hwan[ir]*dr;
	}
	assert(sum>0.0);
	for(int ir=0; ir<nr; ++ir)
	{
		dis_Hwan[ir]/=sum;
	}
	sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_Hwan[ir]*dr;
	}
	cout << "sum is " << sum << endl;

	for(int ir=0; ir<nr; ++ir)
	{
		ofs << dr*(ir+1) << " " << dis_Hwan[ir] << endl;
	}
	ofs.close();
	delete[] dis_Hwan;




	// print out 2D information
	/*
	ofstream ofs2D("dis_wan_2D.dat");

	sum=0.0;
    for(int iy=0; iy<ny; ++iy)
    {
        for(int ix=0; ix<nx; ++ix)
        {
			sum += coord_xy[ix][iy] * INPUT.dx * INPUT.dy;
		}
	}
		

	if(sum>0.0)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			for(int ix=0; ix<nx; ++ix)
			{
				ofs2D << coord_xy[ix][iy]/sum << " ";
			}
			ofs2D << endl;
		}
	}

    for(int ix=0; ix<nx; ++ix)
    {
        delete[] coord_xy[ix];
    }
    delete[] coord_xy;
	
	ofs2D.close();
	*/

	// clean up 
	delete[] dis_total_dipole;	

	if(init_previous)
	{
		delete[] previous_d;
		delete[] previous_dx;
		delete[] previous_dy;
		delete[] previous_dz;
		delete[] previous2_d;
		delete[] previous2_dx;
		delete[] previous2_dy;
		delete[] previous2_dz;
	}
	ofs_dipole.close();
	ofs_vdipole.close();
	ofs_tdipole.close();


	return;
}


//----------------------------------------------------------
// read Wannier functions and calculate the dipole moments
//----------------------------------------------------------
void Wannier::read_wan(Cell &cel, double* dis_Hwan, ofstream &ofs_dipole, ofstream &ofs_vdipole, ofstream &ofs_tdipole)
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
	if(INPUT.ntype==3)
	{ 
		cout << "ito = " << ito << endl;
		cout << "ith = " << ith << endl;
		cout << "itc = " << itc << endl;
		assert(itc>=0); 
	}

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	Water *water = new Water[cel.atom[ito].na];
	Water::nions = 0;
	HBs::setup_water(cel, water);


	if(!init_previous)
	{
		previous_t = -1.0;
		previous2_t = -1.0; 
		previous_d = new double[cel.atom[ito].na]();
		previous_dx = new double[cel.atom[ito].na]();	
		previous_dy = new double[cel.atom[ito].na]();	
		previous_dz = new double[cel.atom[ito].na]();	
		previous2_d = new double[cel.atom[ito].na]();
		previous2_dx = new double[cel.atom[ito].na]();	
		previous2_dy = new double[cel.atom[ito].na]();	
		previous2_dz = new double[cel.atom[ito].na]();	
		init_previous=true;
	}

	bool single_ion = false;
	int ion_index = -1;
	if(INPUT.func == 2)
	{
		assert(INPUT.system == "hydroxide" or INPUT.system == "hydronium");
		if (Water::nions == 1)
		{
			single_ion = true;
			for (int ia=0; ia<cel.atom[ito].na; ia++)
			{
				if((water[ia].nH == 1 and INPUT.system == "hydroxide") or (water[ia].nH == 3 and INPUT.system == "hydronium"))
				{
					ion_index = ia;
				}
			}
		}
	}// renxi 20211120
	
	ofs_dipole << cel.snapshot_index << " " << cel.snapshot_time << endl;
	// change ps to 1.0e-17 means *100000
	
	double ddt = 0.0;
	if(previous2_t>0.0)
	{
		ddt = (cel.snapshot_time - previous2_t)*100000/2.418884326505;
	}
	if(previous_t>0.0)
	{
		previous2_t = previous_t;
	}
	previous_t = cel.snapshot_time;

	// mohan add 2018-02-01
	double total_dipole_x=0.0;
	double total_dipole_y=0.0;
	double total_dipole_z=0.0;
	int iat=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(cel.atom[it].id!=INPUT.ele1) continue;
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			// the first criteria
			if( INPUT.system=="hydronium" and water[ia].nH != 3 and INPUT.func == 1) continue;
			else if( INPUT.system=="hydroxide" and water[ia].nH != 1 and INPUT.func == 1) continue;
			else if( INPUT.system=="water" and water[ia].nH != 2) continue;
			// the second criteria
			if(INPUT.nacc>=0)
			{
				if(water[ia].naccept != INPUT.nacc) continue;
			}
			if(INPUT.ndon>=0)
			{
				if(water[ia].ndonate != INPUT.ndon) continue;
			}
			// the third criteria. renxi added 20211121
			if (INPUT.func == 2 and single_ion== true and ion_index >= 0)
			{
				bool flag = false;
				if (INPUT.system == "hydroxide")
				{
					if (INPUT.func_c == 1)
					{
						for (int iacc=0; iacc<water[ion_index].naccept; iacc++)
						{
							if (ia == water[ion_index].acceptO[iacc])
							{
								flag = true;
							}
						}
					}
					else if (INPUT.func_c == 2)
					{
						for (int idon=0; idon<water[ion_index].ndonate; idon++)
						{
							if (ia == water[ion_index].donateO[idon])
							{
								flag = true;
							}
						}
					}
				}
				if (INPUT.system == "hydronium")
				{
					if (INPUT.func_c == 1)
					{
						for(int idon = 0; idon<water[ion_index].ndonate; idon++)
						{
							if (ia == water[ion_index].donateO[idon])
							{
								flag = true;
							}
						}
					}
					else if (INPUT.func_c == 2)
					{
						for(int iacc = 0; iacc<water[ion_index].naccept; iacc++)
						{
							if (ia == water[ion_index].acceptO[iacc])
							{
								flag = true;
							}
						}
					}
				}
				if (!flag)
				{
					continue;
				}
			}

//			if(INPUT.system == "water" and water[ia].nH !=2 ) continue;

			// mohan add 2018-09-15
			// test
			if(INPUT.system == "water" and water[ia].nH!=2)
			{
				cout << "WARNING!" << endl;
				cout << "water " << ia+1 << " has " << water[ia].nH << " H-bonds" << endl;

				ofs_running << "WARNING!" << endl;
				ofs_running << "water " << ia+1 << " has " << water[ia].nH << " H-bonds" << endl;
				for(int ih=0; ih<water[ia].nH; ++ih)
				{
					int index_test = water[ia].indexH[ih];
					double dx_test=shortest(cel.atom[it].pos[ia].x, cel.atom[ith].pos[index_test].x, norm1);
					double dy_test=shortest(cel.atom[it].pos[ia].y, cel.atom[ith].pos[index_test].y, norm2);
					double dz_test=shortest(cel.atom[it].pos[ia].z, cel.atom[ith].pos[index_test].z, norm3);
					cout << "ih " << ih+1 << " " << sqrt(dx_test*dx_test+dy_test*dy_test+dz_test*dz_test) << endl;
				}
				//exit(0);
			}

			int ih1 = water[ia].indexH[0]; 
		    int ih2 = water[ia].indexH[1]; 

			//-----------------------------------------
			// mohan added 2018-09-15
			// in case there are 3 or more H
			//-----------------------------------------
			if(water[ia].nH>2)
			{
				double* dis_test = new double[ water[ia].nH ];
	
				for(int ih=0; ih<water[ia].nH; ++ih)
				{
					int index_test = water[ia].indexH[ih];
					double dx_test=shortest(cel.atom[it].pos[ia].x, cel.atom[ith].pos[index_test].x, norm1);
     				double dy_test=shortest(cel.atom[it].pos[ia].y, cel.atom[ith].pos[index_test].y, norm2);
					double dz_test=shortest(cel.atom[it].pos[ia].z, cel.atom[ith].pos[index_test].z, norm3);
					dis_test[ih] = sqrt(dx_test*dx_test+dy_test*dy_test+dz_test*dz_test);
				}

				double dis1_test=0.0;
				double dis2_test=0.0;

				if(dis_test[0]<dis_test[1])
				{
					ih1 = water[ia].indexH[0]; dis1_test=dis_test[0];
					ih2 = water[ia].indexH[1]; dis2_test=dis_test[1];
				}
				else
				{
					ih1 = water[ia].indexH[1]; dis1_test=dis_test[1];
					ih2 = water[ia].indexH[0]; dis2_test=dis_test[0];
				}

                for(int ih=0; ih<water[ia].nH; ++ih)
				{
					if(water[ia].indexH[ih] == ih1 or water[ia].indexH[ih] == ih2) continue;

					if(dis_test[ih]>dis1_test and dis_test[ih]>dis2_test) continue;
					else if(dis_test[ih]>dis1_test and dis_test[ih]<dis2_test)
					{
						ih2 = ih; dis2_test=dis_test[ih];
					}
					else if(dis_test[ih]<dis1_test)
					{
						ih2 = ih1; dis2_test=dis1_test;
						ih1 = ih; dis1_test=dis_test[ih];
					}
					else
					{
						cout << "something wrong, please check here." << endl;
						exit(0);
					}
				}

				delete[] dis_test;

				cout << "dis12 = " << dis1_test << " " << dis2_test << endl;
			}

			


			const double dx1 = shortest(cel.atom[it].pos[ia].x, cel.atom[ith].pos[ih1].x, norm1);
			const double dy1 = shortest(cel.atom[it].pos[ia].y, cel.atom[ith].pos[ih1].y, norm2);
			const double dz1 = shortest(cel.atom[it].pos[ia].z, cel.atom[ith].pos[ih1].z, norm3);
			const double dx2 = shortest(cel.atom[it].pos[ia].x, cel.atom[ith].pos[ih2].x, norm1);
			const double dy2 = shortest(cel.atom[it].pos[ia].y, cel.atom[ith].pos[ih2].y, norm2);
			const double dz2 = shortest(cel.atom[it].pos[ia].z, cel.atom[ith].pos[ih2].z, norm3);
			water[ia].dipole[0] = dx1 + dx2;
			water[ia].dipole[1] = dy1 + dy2;
			water[ia].dipole[2] = dz1 + dz2;

//			ofs_running << "dx1,dy1,dz1 = " << dx1 << " " << dy1 << " " << dz1 << endl;
//			ofs_running << "dx2,dy2,dz2 = " << dx2 << " " << dy2 << " " << dz2 << endl;

			// check by mohan
			if( (isnan(water[ia].dipole[0]) or isnan(water[ia].dipole[1]) or isnan(water[ia].dipole[2])) and INPUT.system=="water" )// fixed by renxi 20200414
			{
				cout << "ih1=" << ih1 << " ih2=" << ih2 << endl;
				cout << dx1 << " " << dy1 << " " << dz1 << " " << dx2 << " " << dy2 << " " << dz2 << endl;
				cout << cel.atom[it].pos[ia].y << " " << cel.atom[ith].pos[ih2].y << " " << norm2 << endl;
				cout << water[ia].dipole[0] << " " << water[ia].dipole[1] << " " << water[ia].dipole[2] << endl;
				cout << "warning! dipole is nan!" << endl;
				exit(0);
			}

			for(int ib=0; ib<INPUT.nbands; ++ib)
			{
				// distance between Wannier centers and oxygen atoms
				double dx=shortest(cel.atom[it].pos[ia].x, cel.wan_centers[ib].x, norm1);	
				double dy=shortest(cel.atom[it].pos[ia].y, cel.wan_centers[ib].y, norm2);	
				double dz=shortest(cel.atom[it].pos[ia].z, cel.wan_centers[ib].z, norm3);
				double dis = sqrt(dx*dx+dy*dy+dz*dz);
//				cout << "dis = " << dis << endl;
				assert(dis>=0.0); //added by renxi 20200414
				if(dis < 0.8) // 0.8 Angstroms is a safe cutoff
				{
					//ofs_running << "dx,dy,dz = " << dx << " " << dy << " " << dz << endl;
	//				ofs_running << "wan_centers[" << ib << "]=" 
	//				<< cel.wan_centers[ib].x << " " << cel.wan_centers[ib].y << " " << cel.wan_centers[ib].z << endl;
					water[ia].dipole[0] -= 2.0*dx;
					water[ia].dipole[1] -= 2.0*dy;
					water[ia].dipole[2] -= 2.0*dz;
				}

				// for 1D
				if(dis<INPUT.rcut) 
				{
					int index=dis/this->dr;
					assert(index<this->nr);
					dis_Hwan[index]++;
//					cout << "Wannier center index " << ib+1 << endl;
//					cout << "atom " << iat+1 << " dis " << dis << " Angstroms" << endl;
				}

				// for 2D plot
				/*
				assert(cel.volume>0);
				double density = cel.atom[ito].na*18*1.6605/cel.volume;
				int indexX = (density-INPUT.x0)/INPUT.dx;
				int indexY = (dis - INPUT.y0)/INPUT.dy;
				
				if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0)
				{
					coord_xy[indexX][indexY]+=1.0;
				}
				else
				{
				//	cout << "index=" << indexX << " " << indexY << endl;
				}
				*/
			}// ib

			water[ia].dipole[0] /= -0.20819434;
			water[ia].dipole[1] /= -0.20819434;
			water[ia].dipole[2] /= -0.20819434;


			// change from eA to Deybe
			water[ia].dipole_sum = sqrt(water[ia].dipole[0]*water[ia].dipole[0]
			+water[ia].dipole[1]*water[ia].dipole[1]
			+water[ia].dipole[2]*water[ia].dipole[2]); // 1Debye=0.20819434eA
			ofs_running << "water_dipole " << ia+1 << " " << cel.snapshot_time << " " << water[ia].dipole_sum << endl;

			if(water[ia].dipole_sum>5.0 and INPUT.system=="water")
			{
				cout << "dipole moment wrong! " << water[ia].dipole_sum << endl;
//				exit(0);//commented by renxi 20200414
			}

			// sum of dipole
			total_dipole_x += water[ia].dipole[0];
			total_dipole_y += water[ia].dipole[1];
			total_dipole_z += water[ia].dipole[2];


			if( isnan(water[ia].dipole_sum) and INPUT.system == "water" )
			{
				cout << "warning! the dipole sum is " << water[ia].dipole_sum 
				<< water[ia].dipole[0] << " "
				<< water[ia].dipole[1] << " "
				<< water[ia].dipole[2] << endl;
				exit(0);
			}

			ofs_dipole << water[ia].dipole[0] << " " 
						<< water[ia].dipole[1] << " "
						<< water[ia].dipole[2] << " " 
						<< water[ia].dipole_sum << endl;

			if(water[ia].dipole_sum<INPUT.rcut1) 
			{
				avg_dipole += water[ia].dipole_sum;
				count_dipole++; 
				int index=water[ia].dipole_sum/INPUT.dz;
				assert(index<this->nd);
				dipole_dis[index]++;
			}
			
			// calculating the velocity
			if(ddt>0.0)
			{
				double vdx = (water[ia].dipole[0]-previous2_dx[ia])/ddt;
				double vdy = (water[ia].dipole[1]-previous2_dy[ia])/ddt;
				double vdz = (water[ia].dipole[2]-previous2_dz[ia])/ddt;
				ofs_vdipole << vdx << " " << vdy << " " << vdz << endl;
			}

			if(previous2_t>0.0)
			{
				previous2_dx[ia] = previous_dx[ia]; 
				previous2_dy[ia] = previous_dy[ia];
				previous2_dz[ia] = previous_dz[ia];
				previous2_d[ia] = previous_d[ia];
			}

			if(previous_t>0.0)
			{
				previous_dx[ia] = water[ia].dipole[0]; 
				previous_dy[ia] = water[ia].dipole[1]; 
				previous_dz[ia] = water[ia].dipole[2]; 
				previous_d[ia] = water[ia].dipole_sum;	
			}


			++iat;
		}// ia
	}// it

	// mohan added 2018-02-01
	// total dipole for this snapshot
	double total_dipole = total_dipole_x*total_dipole_x 
		+ total_dipole_y * total_dipole_y
		+ total_dipole_z * total_dipole_z;

	total_dipole /= INPUT.natom1;
	total_dipole = sqrt(total_dipole);

	if(total_dipole<rcut_total_dipole)
	{
		int index = total_dipole/dr_total_dipole;
		dis_total_dipole[index]++;
	}


//	ofs_tdipole << total_dipole << endl; 

	delete[] water;

	return;
}
