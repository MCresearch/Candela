#include "cellFile.h"
#include "input.h"
#include "ili.h"
#include "math.h"

void ILI::Routine()
{
	TITLE("ILI","Routine");
	
	ofs_running << "Compute Instantaneous Liquid Interface (ILI)" << endl;
	ofs_running << "The size of ILI is " << INPUT.nx << " * " << INPUT.ny << endl;

	assert(INPUT.nx>0);
	assert(INPUT.ny>0);

	// setup parameters
//	double vol = INPUT.celldm1*INPUT.celldm2*INPUT.celldm3;
//	this->ref_rho = INPUT.natom1/2.0/vol; // assume the first type of atom is O
	this->ref_rho = INPUT.ref_rho;
	
	ofs_running << "reference rho is " << ref_rho << endl;


	// setup interface
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

	//cout << "geo1=" << INPUT.geo_1 << endl;
	//cout << "geo2=" << INPUT.geo_2 << endl;

	// output data
	ofstream ofs("ILI.dat");

	cout << "------------------------------------------------------------------------------" << endl;
	cout << setw(12) << "Snapshot" << setw(12) << "AvgIter" << setw(15) << "Error(%)" << setw(12) << "NXY" << endl;
	cout << "------------------------------------------------------------------------------" << endl;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{

		// cel : input geometry file
		CellFile cel;

		//ofs_running << "geometry " << igeo%INPUT.geo_interval << endl;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		//cout << igeo << " use:" << cel.read_and_used << endl;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		ofs_running << setw(12) << igeo;
		cout << setw(12) << igeo << " ";

		//get_ili( cel );
		get_ili2( cel );

		// mpi parallel codes
#ifdef __MPI
		int nxy = INPUT.nx*INPUT.ny;
		double* tmp = new double[nxy]();
		double* tmp2 = new double[nxy]();
		// for interface
		for(int ix=0; ix<INPUT.nx; ++ix)
			for(int iy=0; iy<INPUT.ny; ++iy)
				tmp[ix*INPUT.ny+iy]=interface[ix][iy];
		MPI_Allreduce(tmp, tmp2, nxy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		for(int ix=0; ix<INPUT.nx; ++ix)
			for(int iy=0; iy<INPUT.ny; ++iy)
				interface[ix][iy]=tmp2[ix*INPUT.ny+iy];
		// for gradient
		for(int j=0; j<3; ++j)
		{
			for(int ix=0; ix<INPUT.nx; ++ix)
				for(int iy=0; iy<INPUT.ny; ++iy)
					tmp[ix*INPUT.ny+iy]=gradient[ix][iy][j];
			MPI_Allreduce(tmp, tmp2, nxy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			for(int ix=0; ix<INPUT.nx; ++ix)
				for(int iy=0; iy<INPUT.ny; ++iy)
					gradient[ix][iy][j]=tmp2[ix*INPUT.ny+iy];
		}
		delete[] tmp;
		delete[] tmp2;
#endif

		// print out the information
		if(RANK==0)
		{
			ofs << "ILI" << igeo << endl;
			ofs << INPUT.nx << " " << INPUT.ny << endl;
			for(int ix=0; ix<INPUT.nx; ++ix)
			{
				for(int iy=0; iy<INPUT.ny; ++iy)
				{
					ofs << interface[ix][iy] << " ";
				}
				ofs << endl;
			}
			for(int j=0; j<3; ++j)
			{
				ofs << "GRAD" << j+1 << endl;
				ofs << INPUT.nx << " " << INPUT.ny << endl;
				for(int ix=0; ix<INPUT.nx; ++ix)
				{
					for(int iy=0; iy<INPUT.ny; ++iy)
					{
						ofs << gradient[ix][iy][j] << " ";
					}
					ofs << endl;
				}
			}
		}
	}	
	ofs.close();

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

	return;
}


// the new algorithm to find the instantanous liquid interface,
// mohan wrote on 2017-08-02
void ILI::get_ili2(const Cell &cel)
{
	TITLE("ILI","get_ili2");

	assert( INPUT.zeta > 0.0 );

	// prepare data
	int ntry = INPUT.ntry;
	double z0 = INPUT.z0;
	double dz = INPUT.dz;

	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			interface[ix][iy]=0.0;
			for(int j=0; j<3; ++j)
			{
				this->gradient[ix][iy][j]=0.0;
			}
		}
	}
	
	// prepare data 2
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double dx = norm1/INPUT.nx;
	double dy = norm2/INPUT.ny;

	int ixy=0;
	int sxy=0;
	int avg_iter=0;
	double avg_error=0.0;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
#ifdef __MPI
			++ixy;
			if( ixy%NPROC!=RANK ) continue;
#endif

			double outside_z=z0;
			double z_boundary=0.0;
			double move_z=0.0;
			double drho=0.0;
			int iter=0;
			double z_now=0.0;

calculate_rho:

			// we approach the surface from the vacuum,
			// so drho=0.0 means the first step
			// while drho<0.0 means we are still in the vacuum
			// and z_boundary=0.0 means we are not searching the ILI
			// beyond z_boundary.
			if( drho==0.0 or (drho<0.0 and z_boundary==0.0) )
			{
				move_z += dz;
				z_now = outside_z + move_z;
			}
			// We are still approaching the surface but we already
			// know the boundary.
			else if(drho<0.0 and z_boundary>0.0)
			{
				z_now = (z_now+z_boundary)*0.5;
			}
			// We are crossing the boundary because drho>0.0
			// so we set new boundary (should be closer to the point we want)
			// and recalculate the z_now.
			else if(drho>0.0) 
			{
				z_boundary = z_now;
				double z_previous = z_now - move_z; 
				z_now = (z_boundary + z_previous)*0.5;
			}

			// false means do not calculate the gradient
			double rho_now = cal_gauss(cel, ix*dx, iy*dy, z_now, false, gradient[ix][iy]);
			//cout << setw(20) << z_now << setw(20) << rho_now << endl;

			if(rho_now < (1.00+INPUT.within)*ref_rho and rho_now > (1.00-INPUT.within)*ref_rho)
			{
				// true means calculate the gradient
				cal_gauss(cel, ix*dx, iy*dy, z_now, true, gradient[ix][iy]);
				interface[ix][iy] = z_now;
				avg_iter+=iter;
				if(ref_rho>0) avg_error+=abs((rho_now-ref_rho)/ref_rho);
				++sxy;
			}
			else if(iter>=INPUT.maxiter)
			{
				cout << "Target rho not Found! Gradient set to 0. Please check later." << endl;
				ofs_running << "Target rho not Found! Gradient set to 0. Please check later." << endl;
			}
			else
			{
				drho = rho_now-ref_rho;
				++iter;
				goto calculate_rho;
			}
		}//end iy
	}//end ix
	
	// print out statistics
	if(sxy>0)
	{
		double iter_per_snapshot = (double)avg_iter/sxy;
		double error_per_snapshot = avg_error/(double)sxy; 
		ofs_running << " average iteration number is " << iter_per_snapshot << endl;
		ofs_running << " average error is (%) " << error_per_snapshot*100 << endl;
		ofs_running << " successful operation number is " << sxy << endl;
		cout << setw(12) << iter_per_snapshot << setw(15) << error_per_snapshot*100 << setw(12) << sxy << endl;
	}


	return;
}

double ILI::cal_gauss(const Cell &cel, const double &x_in, const double &y_in, const double &z_in, const bool &cal_grad,
double* grad)
{
	// setup cell parameters
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();

	// setup gaussian parameters
	assert( INPUT.zeta > 0.0 );
	double zeta = INPUT.zeta; 
	double zeta2 = zeta*zeta;
    double zeta22 = 2*zeta*zeta;
	double d = (double)INPUT.d;
    double fac = pow(3.1415926*zeta22, -d/2);

	double rho=0.0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(INPUT.ele1!="none" and cel.atom[it].id!=INPUT.ele1) {continue;}
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			double posx = cel.atom[it].pos[ia].x;
			double posy = cel.atom[it].pos[ia].y;
			double posz = cel.atom[it].pos[ia].z;
			while( posx >= norm1 ) posx -= norm1;	
			while( posx < 0 ) posx += norm1;
			while( posy >= norm2 ) posy -= norm2;	
			while( posy < 0 ) posy += norm2;
			double drx = shortest(x_in,posx,INPUT.celldm1); 
			double dry = shortest(y_in,posy,INPUT.celldm2);
			double drz = 1000.0;
			if(INPUT.system=="water")
			{
				while( posz >= norm3 ) posz -= norm3;	
				while( posz < 0.0 ) posz += norm3;	
			}
			else 
			{
				while( posz > INPUT.upper_z ) posz -= norm3;	
			}
			drz = z_in - posz;

/*
			double r2 = drx*drx + dry*dry + drz*drz;
			double gaussian = fac*exp(-r2/zeta22);

			if(cal_grad)
			{
				grad[0]+=gaussian*(-drx)/zeta2;
				grad[1]+=gaussian*(-dry)/zeta2;
				grad[2]+=gaussian*(-drz)/zeta2;
			}

			rho+=gaussian;	
*/


			for(int ipx=-1; ipx<=1; ++ipx)
			{
				for(int ipy=-1; ipy<=1; ++ipy)
				{
					double new_drx = drx + ipx*INPUT.celldm1;
					double new_dry = dry + ipy*INPUT.celldm2;
					double new_drz = drz;

					//new_drx = abs(new_drx);
					//new_dry = abs(new_dry);
					//new_drz = abs(new_drz);

					double r2 = new_drx*new_drx + new_dry*new_dry + new_drz*new_drz;
					double gaussian = fac*exp(-r2/zeta22);

					if(cal_grad)
					{
						grad[0]+=gaussian*(-new_drx)/zeta2;
						grad[1]+=gaussian*(-new_dry)/zeta2;
						grad[2]+=gaussian*(-new_drz)/zeta2;
					}

					rho+=gaussian;	
				}//end ipy
			}//end ipx
		}
	}

	return rho;
}

void ILI::get_ili(const Cell &cel)
{
	TITLE("ILI","get_ili");

	assert( INPUT.zeta > 0.0 );

	// prepare data
	int ntry = INPUT.ntry;
	double z0 = INPUT.z0;
	double dz = INPUT.dz;
	double* try_z = new double[ntry]; // y
	double* rho_z = new double[ntry]; // x

	for(int iz=0; iz<ntry; ++iz)
	{
		try_z[iz]=z0+dz*iz;
	}

	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
			interface[ix][iy]=0.0;
			for(int j=0; j<3; ++j)
			{
				this->gradient[ix][iy][j]=0.0;
			}
		}
	}
	
	// prepare data 2
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double dx = norm1/INPUT.nx;
	double dy = norm2/INPUT.ny;

	int ixy=0;
	double lowest_z=1000.0; // mohan test
	double highest_z=-1000.0; // mohan test
	int x0=-1;
	int y0=-1;
	int x1=-1;
	int y1=-1;
	for(int ix=0; ix<INPUT.nx; ++ix)
	{
		for(int iy=0; iy<INPUT.ny; ++iy)
		{
#ifdef __MPI
			++ixy;
			if( ixy%NPROC!=RANK ) continue;
#endif
			//ofs_running << setw(10) << ix << setw(10) << iy << endl;
			for(int iz=0; iz<ntry; ++iz)
			{
				// when preparing for interpolation, one does not need to calculate gradient 
				rho_z[iz] = cal_gauss(cel, ix*dx, iy*dy, try_z[iz], false, gradient[ix][iy]);
			}

			// use x and y arrays to do interpolation!
			double target_z = Polynomial_Interpolation_xy(rho_z, try_z, ntry, ref_rho);
			if(target_z<0 or target_z>=INPUT.celldm3)
			{
				ofs_running << "something is wrong." << endl;
				ofs_running << setw(20) << "computed_rho" << setw(20) << "z" << endl; 
				for(int iz=0; iz<ntry; ++iz)
				{
					ofs_running << setw(20) << rho_z[iz] << setw(20) << try_z[iz] <<endl;
				}
				ofs_running << "ref_rho=" << ref_rho << endl;
				ofs_running << "target_z=" << target_z << endl;
				target_z=0.0;
			}
			interface[ix][iy] = target_z;
			if(lowest_z>target_z)
			{
				lowest_z=target_z;
				x0 = ix;	
				y0 = iy;
			}
			if(highest_z,target_z)
			{
				highest_z=target_z;
				x1 = ix;	
				y1 = iy;
			}

			// gradient is calculated here
			double actual_rho = cal_gauss(cel, ix*dx, iy*dy, target_z, true, gradient[ix][iy]);
			if(actual_rho > 1.03*ref_rho or actual_rho < 0.97*ref_rho)
			{
//				for(int iz=0; iz<ntry; ++iz)
//				{
//					ofs_running << setw(15) << try_z[iz] << setw(15) << rho_z[iz] << endl;
//				}
				ofs_running << "the interpolated z is " << target_z << endl;
				ofs_running << "warning! the rho from interpolated z is " << actual_rho << ", however the reference rho is " << ref_rho << endl;
			}
		}
	}


	cout << setw(15) << "highest_z" << setw(15) << highest_z 
		<< setw(5) << "x1" << setw(5) << x1 
		<< setw(5) << "y1" << setw(5) << y1;

	cout << setw(15) << "lowest_z" <<  setw(15) << lowest_z 
		<< setw(5) << "x0" << setw(5) << x0 
		<< setw(5) << "y0" << setw(5) << y0 << endl;


	delete[] try_z;
	delete[] rho_z;
	return;
}
